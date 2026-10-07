# PyGEL cheatsheet

Week 01 uses `pygel3d.hmesh` (halfedge `Manifold`). Indices are ints. Circulation is **CCW**.

```python
from pygel3d import hmesh, jupyter_display as jd
import numpy as np

jd.set_export_mode(True)
m = hmesh.load("bunnygtest.obj")   # OBJ/OFF/PLY/X3D; None on fail
hmesh.save("out.obj", m)
jd.display(m)
m2 = hmesh.Manifold(m)             # copy
```

## Counts (holes in arrays)

`no_allocated_*` is **allocation size**, not live count. After deletes, call `m.cleanup()` (invalidates extra arrays).

```python
m.no_allocated_vertices()
m.no_allocated_faces()
m.no_allocated_halfedges()
m.vertices()   # iterable of live vids
m.faces()      # iterable of live fids
m.halfedges()
```

## Positions

`pos = m.positions()` is a **view**, not a copy. Edits write through.

```python
pos = m.positions()          # (n,3)
pos[vid]                     # vertex vid
pos[:] = new_pos             # write all
```

## Circulate (the important API)

```python
# around vertex vid
m.circulate_vertex(vid, mode="v")  # neighbour vertices (one-ring)
m.circulate_vertex(vid, mode="f")  # incident faces
m.circulate_vertex(vid, mode="h")  # outgoing halfedges

# around face fid
m.circulate_face(fid, mode="v")    # vertices of face
m.circulate_face(fid, mode="f")    # adjacent faces (edge-neighbours)
m.circulate_face(fid, mode="h")    # halfedges of face
```

Always `list(...)` if you need random access / length.

```python
nbrs = list(m.circulate_vertex(vid, mode="v"))
fverts = list(m.circulate_face(fid, mode="v"))
```

## Halfedges

Each edge = two opposite halfedges.

| call                     | meaning                      |
| ------------------------ | ---------------------------- |
| `m.next_halfedge(h)`     | next around face             |
| `m.prev_halfedge(h)`     | prev around face             |
| `m.opposite_halfedge(h)` | other side of edge           |
| `m.incident_face(h)`     | face of `h`                  |
| `m.incident_vertex(h)`   | vertex **pointed to** by `h` |

Vertex walk: `h → opposite → next`. Face walk: follow `next`.

```python
m.valency(vid)                 # |one-ring|
m.connected(v0, v1)
m.face_normal(fid)
m.vertex_normal(vid)           # angle-weighted
m.no_edges(fid)                # edges of that face
```

## Build faces (dual / new mesh)

```python
m2 = hmesh.Manifold()
fid = m2.add_face([p0, p1, p2, ...])   # list of 3D points → new face
# or
m2 = hmesh.Manifold.from_triangles(verts, tris)
```

`add_face` **creates new vertices from coordinates**. Shared vertices only if you reuse the same point objects carefully — for dual, usually: one new vertex per old face (centroid), then `add_face` with those centroids in CCW order from `circulate_vertex(vid, mode="f")`.

## Built-ins you can use / check against

```python
hmesh.volume(m)            # signed; mesh must be closed
hmesh.area(m)
hmesh.closed(m)
hmesh.valid(m)
hmesh.flip_orientation(m)
hmesh.triangulate(m)
hmesh.laplacian_smooth(m, w, iters)
```

## Week 01 patterns

**Laplacian smooth** — average one-ring, write after a full pass (Jacobi):

```python
pos = m.positions()
new_pos = np.zeros_like(pos)
for vid in m.vertices():
    nbrs = list(m.circulate_vertex(vid, mode="v"))
    new_pos[vid] = pos[nbrs].mean(axis=0)
pos[:] = new_pos
```

**Volume** — signed tets with origin. For each face, vertices \(p_0,p_1,\ldots\) in CCW:

\[
V=\sum_f \frac16\, p_0\cdot(p_1\times p_2)
\]

(triangle faces). Sign follows orientation; take `abs` if you only want size.

```python
pos = m.positions()
vol = 0.0
for fid in m.faces():
    vs = list(m.circulate_face(fid, mode="v"))
    p0, p1, p2 = pos[vs[0]], pos[vs[1]], pos[vs[2]]
    vol += np.dot(p0, np.cross(p1, p2)) / 6.0
```

**Dual** — new vertex ↔ old face; new face ↔ old vertex.

1. For each old `fid`, store centroid of `pos[circulate_face(fid,"v")]`.
2. For each old `vid`, `add_face` the centroids of `circulate_vertex(vid, mode="f")` (CCW).
