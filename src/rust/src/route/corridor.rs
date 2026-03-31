use std::collections::BinaryHeap;
use std::time::Instant;

// ---------------------------------------------------------------------------
// CSR (Compressed Sparse Row) graph
// ---------------------------------------------------------------------------

struct CsrGraph {
    offsets: Vec<usize>,
    targets: Vec<usize>,
    weights: Vec<f64>,
}

impl CsrGraph {
    fn n_nodes(&self) -> usize {
        self.offsets.len() - 1
    }

    fn n_edges(&self) -> usize {
        self.targets.len()
    }
}

// ---------------------------------------------------------------------------
// BinaryHeap min-heap state
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct State {
    cost: f64,
    node: usize,
}

impl PartialEq for State {
    fn eq(&self, other: &Self) -> bool {
        self.cost == other.cost && self.node == other.node
    }
}
impl Eq for State {}

impl Ord for State {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Reverse: smaller cost = higher priority; break ties by node index
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| other.node.cmp(&self.node))
    }
}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

// ---------------------------------------------------------------------------
// Neighbor offset tables
// ---------------------------------------------------------------------------

struct NeighborOffset {
    dr: i32,
    dc: i32,
    dist_fn: DistKind,
}

enum DistKind {
    CellWidth,
    CellHeight,
    Diagonal,
    Knight(i32, i32), // (|dr|, |dc|) for distance computation
}

fn neighbor_offsets(neighbours: i32) -> Vec<NeighborOffset> {
    let mut offsets = vec![
        // 4-connectivity (cardinal)
        NeighborOffset { dr: 0, dc: -1, dist_fn: DistKind::CellWidth },
        NeighborOffset { dr: 0, dc: 1, dist_fn: DistKind::CellWidth },
        NeighborOffset { dr: -1, dc: 0, dist_fn: DistKind::CellHeight },
        NeighborOffset { dr: 1, dc: 0, dist_fn: DistKind::CellHeight },
    ];

    if neighbours >= 8 {
        // Diagonals
        offsets.push(NeighborOffset { dr: -1, dc: -1, dist_fn: DistKind::Diagonal });
        offsets.push(NeighborOffset { dr: -1, dc: 1, dist_fn: DistKind::Diagonal });
        offsets.push(NeighborOffset { dr: 1, dc: -1, dist_fn: DistKind::Diagonal });
        offsets.push(NeighborOffset { dr: 1, dc: 1, dist_fn: DistKind::Diagonal });
    }

    if neighbours >= 16 {
        // Knight's moves
        offsets.push(NeighborOffset { dr: -2, dc: -1, dist_fn: DistKind::Knight(2, 1) });
        offsets.push(NeighborOffset { dr: -2, dc: 1, dist_fn: DistKind::Knight(2, 1) });
        offsets.push(NeighborOffset { dr: 2, dc: -1, dist_fn: DistKind::Knight(2, 1) });
        offsets.push(NeighborOffset { dr: 2, dc: 1, dist_fn: DistKind::Knight(2, 1) });
        offsets.push(NeighborOffset { dr: -1, dc: -2, dist_fn: DistKind::Knight(1, 2) });
        offsets.push(NeighborOffset { dr: -1, dc: 2, dist_fn: DistKind::Knight(1, 2) });
        offsets.push(NeighborOffset { dr: 1, dc: -2, dist_fn: DistKind::Knight(1, 2) });
        offsets.push(NeighborOffset { dr: 1, dc: 2, dist_fn: DistKind::Knight(1, 2) });
    }

    offsets
}

fn offset_distance(kind: &DistKind, cell_w: f64, cell_h: f64, diag: f64) -> f64 {
    match kind {
        DistKind::CellWidth => cell_w,
        DistKind::CellHeight => cell_h,
        DistKind::Diagonal => diag,
        DistKind::Knight(dr, dc) => {
            let dy = (*dr as f64) * cell_h;
            let dx = (*dc as f64) * cell_w;
            (dy * dy + dx * dx).sqrt()
        }
    }
}

// ---------------------------------------------------------------------------
// Graph construction (two-pass CSR)
// ---------------------------------------------------------------------------

/// Build a CSR graph from a raster grid. Returns (graph, min_friction).
fn build_grid_graph(
    values: &[f64],
    n_rows: usize,
    n_cols: usize,
    cell_w: f64,
    cell_h: f64,
    neighbours: i32,
) -> (CsrGraph, f64) {
    let n_cells = n_rows * n_cols;
    let offsets_table = neighbor_offsets(neighbours);
    let diag = (cell_w * cell_w + cell_h * cell_h).sqrt();

    // Precompute distances for each offset
    let distances: Vec<f64> = offsets_table
        .iter()
        .map(|o| offset_distance(&o.dist_fn, cell_w, cell_h, diag))
        .collect();

    // Pass 1: count edges per node
    let mut degree = vec![0u32; n_cells];
    let mut min_friction = f64::INFINITY;

    for r in 0..n_rows {
        for c in 0..n_cols {
            let idx = r * n_cols + c;
            let v = values[idx];
            if v.is_nan() || !v.is_finite() {
                continue;
            }
            if v < min_friction {
                min_friction = v;
            }
            for off in &offsets_table {
                let nr = r as i32 + off.dr;
                let nc = c as i32 + off.dc;
                if nr < 0 || nr >= n_rows as i32 || nc < 0 || nc >= n_cols as i32 {
                    continue;
                }
                let nv = values[nr as usize * n_cols + nc as usize];
                if nv.is_nan() || !nv.is_finite() {
                    continue;
                }
                degree[idx] += 1;
            }
        }
    }

    // Build offsets via prefix sum
    let mut csr_offsets = vec![0usize; n_cells + 1];
    for i in 0..n_cells {
        csr_offsets[i + 1] = csr_offsets[i] + degree[i] as usize;
    }
    let total_edges = csr_offsets[n_cells];

    let mut targets = vec![0usize; total_edges];
    let mut weights = vec![0.0f64; total_edges];

    // Pass 2: fill edges
    let mut pos = vec![0u32; n_cells]; // write cursor per node

    for r in 0..n_rows {
        for c in 0..n_cols {
            let idx = r * n_cols + c;
            let v = values[idx];
            if v.is_nan() || !v.is_finite() {
                continue;
            }
            for (oi, off) in offsets_table.iter().enumerate() {
                let nr = r as i32 + off.dr;
                let nc = c as i32 + off.dc;
                if nr < 0 || nr >= n_rows as i32 || nc < 0 || nc >= n_cols as i32 {
                    continue;
                }
                let nidx = nr as usize * n_cols + nc as usize;
                let nv = values[nidx];
                if nv.is_nan() || !nv.is_finite() {
                    continue;
                }
                let w = ((v + nv) / 2.0) * distances[oi];
                let slot = csr_offsets[idx] + pos[idx] as usize;
                targets[slot] = nidx;
                weights[slot] = w;
                pos[idx] += 1;
            }
        }
    }

    let graph = CsrGraph {
        offsets: csr_offsets,
        targets,
        weights,
    };

    (graph, min_friction)
}

// ---------------------------------------------------------------------------
// Path reconstruction
// ---------------------------------------------------------------------------

fn reconstruct_path(pred: &[usize], origin: usize, dest: usize) -> Vec<usize> {
    let mut path = Vec::new();
    let mut cur = dest;
    while cur != origin {
        path.push(cur);
        cur = pred[cur];
        if cur == usize::MAX {
            return Vec::new(); // no path
        }
    }
    path.push(origin);
    path.reverse();
    path
}

// ---------------------------------------------------------------------------
// Dijkstra
// ---------------------------------------------------------------------------

fn dijkstra(graph: &CsrGraph, origin: usize, dest: usize) -> Option<(Vec<usize>, f64)> {
    let n = graph.n_nodes();
    let mut dist = vec![f64::INFINITY; n];
    let mut pred = vec![usize::MAX; n];
    let mut heap = BinaryHeap::new();

    dist[origin] = 0.0;
    heap.push(State { cost: 0.0, node: origin });

    while let Some(State { cost, node }) = heap.pop() {
        if node == dest {
            let path = reconstruct_path(&pred, origin, dest);
            return Some((path, cost));
        }
        if cost > dist[node] {
            continue; // stale entry
        }
        let start = graph.offsets[node];
        let end = graph.offsets[node + 1];
        for i in start..end {
            let next = graph.targets[i];
            let w = graph.weights[i];
            let new_cost = cost + w;
            if new_cost < dist[next] {
                dist[next] = new_cost;
                pred[next] = node;
                heap.push(State { cost: new_cost, node: next });
            }
        }
    }

    None // no path
}

// ---------------------------------------------------------------------------
// Bidirectional Dijkstra
// ---------------------------------------------------------------------------

fn bidirectional(graph: &CsrGraph, origin: usize, dest: usize) -> Option<(Vec<usize>, f64)> {
    if origin == dest {
        return Some((vec![origin], 0.0));
    }

    let n = graph.n_nodes();
    let mut dist_f = vec![f64::INFINITY; n];
    let mut dist_b = vec![f64::INFINITY; n];
    let mut pred_f = vec![usize::MAX; n];
    let mut pred_b = vec![usize::MAX; n];
    let mut settled_f = vec![false; n];
    let mut settled_b = vec![false; n];

    let mut heap_f = BinaryHeap::new();
    let mut heap_b = BinaryHeap::new();

    dist_f[origin] = 0.0;
    dist_b[dest] = 0.0;
    heap_f.push(State { cost: 0.0, node: origin });
    heap_b.push(State { cost: 0.0, node: dest });

    let mut mu = f64::INFINITY;
    let mut meeting = usize::MAX;

    loop {
        // Check termination
        let min_f = heap_f.peek().map_or(f64::INFINITY, |s| s.cost);
        let min_b = heap_b.peek().map_or(f64::INFINITY, |s| s.cost);

        if min_f + min_b >= mu {
            break;
        }
        if min_f.is_infinite() && min_b.is_infinite() {
            break; // both exhausted, no path
        }

        // Expand the side with smaller minimum key
        if min_f <= min_b {
            if let Some(State { cost, node }) = heap_f.pop() {
                if cost > dist_f[node] {
                    continue;
                }
                settled_f[node] = true;

                // Check meeting
                if settled_b[node] {
                    let total = dist_f[node] + dist_b[node];
                    if total < mu {
                        mu = total;
                        meeting = node;
                    }
                }

                let start = graph.offsets[node];
                let end = graph.offsets[node + 1];
                for i in start..end {
                    let next = graph.targets[i];
                    let w = graph.weights[i];
                    let new_cost = cost + w;
                    if new_cost < dist_f[next] {
                        dist_f[next] = new_cost;
                        pred_f[next] = node;
                        heap_f.push(State { cost: new_cost, node: next });
                    }
                    // Also check if this neighbor is settled backward
                    if dist_b[next] < f64::INFINITY {
                        let total = new_cost + dist_b[next];
                        if total < mu {
                            mu = total;
                            meeting = next;
                        }
                    }
                }
            }
        } else {
            if let Some(State { cost, node }) = heap_b.pop() {
                if cost > dist_b[node] {
                    continue;
                }
                settled_b[node] = true;

                if settled_f[node] {
                    let total = dist_f[node] + dist_b[node];
                    if total < mu {
                        mu = total;
                        meeting = node;
                    }
                }

                let start = graph.offsets[node];
                let end = graph.offsets[node + 1];
                for i in start..end {
                    let next = graph.targets[i];
                    let w = graph.weights[i];
                    let new_cost = cost + w;
                    if new_cost < dist_b[next] {
                        dist_b[next] = new_cost;
                        pred_b[next] = node;
                        heap_b.push(State { cost: new_cost, node: next });
                    }
                    if dist_f[next] < f64::INFINITY {
                        let total = dist_f[next] + new_cost;
                        if total < mu {
                            mu = total;
                            meeting = next;
                        }
                    }
                }
            }
        }
    }

    if meeting == usize::MAX {
        return None;
    }

    // Reconstruct: forward path origin → meeting, backward path meeting → dest
    let mut path = Vec::new();
    let mut cur = meeting;
    while cur != origin {
        path.push(cur);
        cur = pred_f[cur];
    }
    path.push(origin);
    path.reverse();

    cur = pred_b[meeting];
    while cur != usize::MAX && cur != dest {
        path.push(cur);
        cur = pred_b[cur];
    }
    if meeting != dest {
        path.push(dest);
    }

    Some((path, mu))
}

// ---------------------------------------------------------------------------
// A* search
// ---------------------------------------------------------------------------

#[inline]
fn heuristic(
    cell: usize,
    dest: usize,
    n_cols: usize,
    cell_w: f64,
    cell_h: f64,
    min_friction: f64,
) -> f64 {
    let r1 = (cell / n_cols) as f64;
    let c1 = (cell % n_cols) as f64;
    let r2 = (dest / n_cols) as f64;
    let c2 = (dest % n_cols) as f64;
    let dx = (c1 - c2).abs() * cell_w;
    let dy = (r1 - r2).abs() * cell_h;
    (dx * dx + dy * dy).sqrt() * min_friction
}

fn astar(
    graph: &CsrGraph,
    origin: usize,
    dest: usize,
    n_cols: usize,
    cell_w: f64,
    cell_h: f64,
    min_friction: f64,
) -> Option<(Vec<usize>, f64)> {
    let n = graph.n_nodes();
    let mut g_score = vec![f64::INFINITY; n];
    let mut pred = vec![usize::MAX; n];
    let mut heap = BinaryHeap::new();

    g_score[origin] = 0.0;
    let h = heuristic(origin, dest, n_cols, cell_w, cell_h, min_friction);
    heap.push(State { cost: h, node: origin });

    while let Some(State { cost: _f_cost, node }) = heap.pop() {
        if node == dest {
            let path = reconstruct_path(&pred, origin, dest);
            return Some((path, g_score[dest]));
        }

        let g = g_score[node];
        // Skip stale entries: if the best f for this node was already beaten,
        // skip. We compare g-scores since h is deterministic.
        // A node may appear multiple times; we rely on g_score check below.

        let start = graph.offsets[node];
        let end = graph.offsets[node + 1];
        for i in start..end {
            let next = graph.targets[i];
            let w = graph.weights[i];
            let new_g = g + w;
            if new_g < g_score[next] {
                g_score[next] = new_g;
                pred[next] = node;
                let f = new_g + heuristic(next, dest, n_cols, cell_w, cell_h, min_friction);
                heap.push(State { cost: f, node: next });
            }
        }
    }

    None
}

// ---------------------------------------------------------------------------
// Solver result
// ---------------------------------------------------------------------------

pub struct CorridorResult {
    /// 0-based cell indices along the optimal path (empty if no path)
    pub path_cells: Vec<usize>,
    pub total_cost: f64,
    pub solve_time_ms: f64,
    pub graph_build_time_ms: f64,
    pub n_edges: usize,
}

// ---------------------------------------------------------------------------
// Opaque cached graph (CsrGraph stays private)
// ---------------------------------------------------------------------------

pub struct CorridorGraph {
    graph: CsrGraph,
    min_friction: f64,
    pub n_rows: usize,
    pub n_cols: usize,
    pub cell_w: f64,
    pub cell_h: f64,
    pub neighbours: i32,
    pub build_time_ms: f64,
}

impl CorridorGraph {
    pub fn n_edges(&self) -> usize {
        self.graph.n_edges()
    }

    pub fn memory_bytes_estimate(&self) -> usize {
        (self.graph.offsets.len() + self.graph.targets.len()) * std::mem::size_of::<usize>()
            + self.graph.weights.len() * std::mem::size_of::<f64>()
    }
}

/// Build a cached corridor graph from a flattened raster.
pub fn build_corridor_graph(
    values: &[f64],
    n_rows: usize,
    n_cols: usize,
    cell_w: f64,
    cell_h: f64,
    neighbours: i32,
) -> CorridorGraph {
    let t_build = Instant::now();
    let (graph, min_friction) = build_grid_graph(values, n_rows, n_cols, cell_w, cell_h, neighbours);
    let build_time_ms = t_build.elapsed().as_secs_f64() * 1000.0;

    CorridorGraph {
        graph,
        min_friction,
        n_rows,
        n_cols,
        cell_w,
        cell_h,
        neighbours,
        build_time_ms,
    }
}

/// Route on a pre-built cached graph.
pub fn solve_on_graph(
    cg: &CorridorGraph,
    origin: usize,
    dest: usize,
    method: &str,
) -> CorridorResult {
    let t_solve = Instant::now();
    let result = match method {
        "dijkstra" => dijkstra(&cg.graph, origin, dest),
        "bidirectional" => bidirectional(&cg.graph, origin, dest),
        "astar" => astar(
            &cg.graph, origin, dest,
            cg.n_cols, cg.cell_w, cg.cell_h, cg.min_friction,
        ),
        _ => None,
    };
    let solve_time_ms = t_solve.elapsed().as_secs_f64() * 1000.0;

    match result {
        Some((path, total_cost)) => CorridorResult {
            path_cells: path,
            total_cost,
            solve_time_ms,
            graph_build_time_ms: 0.0,
            n_edges: cg.graph.n_edges(),
        },
        None => CorridorResult {
            path_cells: Vec::new(),
            total_cost: f64::INFINITY,
            solve_time_ms,
            graph_build_time_ms: 0.0,
            n_edges: cg.graph.n_edges(),
        },
    }
}

// ---------------------------------------------------------------------------
// Public solver entry point (builds graph + solves in one call)
// ---------------------------------------------------------------------------

pub fn solve(
    values: &[f64],
    n_rows: usize,
    n_cols: usize,
    cell_w: f64,
    cell_h: f64,
    origin: usize,
    dest: usize,
    neighbours: i32,
    method: &str,
) -> CorridorResult {
    let cg = build_corridor_graph(values, n_rows, n_cols, cell_w, cell_h, neighbours);
    let mut result = solve_on_graph(&cg, origin, dest, method);
    result.graph_build_time_ms = cg.build_time_ms;
    result
}
