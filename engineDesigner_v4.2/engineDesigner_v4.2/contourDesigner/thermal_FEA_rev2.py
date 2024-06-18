import matplotlib.pyplot as plt
import meshpy.triangle as triangle
import numpy as np
import scipy.linalg as la


def round_trip_connect(start, end):
    return [(i, i + 1) for i in range(start, end)] + [(end, start)]

def generate_bdry(corner_points, res):
    bdry = []
    numPts = np.zeros(len(corner_points)) # Number of points each edge will be subdivided into, according to res
    for i in np.arange(len(numPts)):
        # Loop through each pair of adjacent points, including end --> start
        a = np.mod(i + 1, len(numPts)) - 1 # Makes sure you loop back to beginning
        dist = np.linalg.norm(np.array(corner_points[a + 1]) - np.array(corner_points[a])) # Euclidean distance
        numPts[i] = np.ceil(dist / res) + 1
        x = np.linspace(corner_points[a][0], corner_points[a + 1][0], numPts[i].astype('int')) # x indices of subdivided points
        y = np.linspace(corner_points[a][1], corner_points[a + 1][1], numPts[i].astype('int')) # y indices of subdivided points
        for j in range(len(x) - 1):
            bdry.append([x[j], y[j]])

    bdry_faces = round_trip_connect(0, len(bdry) - 1) # Generate edges along boundary
    return (bdry, bdry_faces)

def generate_mesh(corner_points, res):
    (bdry, bdry_faces) = generate_bdry(corner_points, res)
    mesh_info = triangle.MeshInfo() # Generate info for mesh
    mesh_info.set_points(bdry)  # Mesh boundary points
    mesh_info.set_facets(bdry_faces) # Mesh boundary edges

    mesh = triangle.build(mesh_info, generate_faces=True, volume_constraints=True,
                          max_volume=res ** 2) # Build mesh - note max triangle area set to res^2 (e.g. a square of side length res)

    mesh_points = np.array(mesh.points)
    mesh_faces = np.array(mesh.faces)
    mesh_tris = np.array(mesh.elements)


    plt.figure()
    plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris)
    plt.xlabel("x")
    plt.ylabel("y")


    return (mesh, mesh_points, mesh_faces, mesh_tris, bdry, bdry_faces)

class MeshPoint:
    def __init__(self, mesh, idx, coords, adj, rel_coords, areas, tris, eq, bdry_bool):
        self.mesh = mesh
        self.mesh_points = np.array(mesh.points)
        self.idx = idx
        self.coords = coords
        self.adj = adj
        self.rel_coords = rel_coords
        self.areas = areas
        self.tris = tris
        self.eq = eq
        self.bdry_bool = bdry_bool
        return

    def heat_balance(self, k, bdry_cond):
        pass
        return


def fea_solver(corner_points, res, boundary_conds, k):

    #Generate mesh from boundary
    (mesh, mesh_points, mesh_faces, mesh_tris, bdry, bdry_facets) = generate_mesh(corner_points, res)

    adj_arr = [[] for i in range(len(mesh_points))]
    for i in range(len(mesh_faces)):
        # For each pair of adjacent points
        # List idx point 1 as a neighbor of point 2, and vice versa
        adj_arr[mesh_faces[i][0]].append(mesh_faces[i][1])
        adj_arr[mesh_faces[i][1]].append(mesh_faces[i][0])

    tri_arr = [[] for i in range(len(mesh_points))]
    for i in range(len(mesh_tris)):
        # For each triangle
        # Add index of triangle in mesh_tris to corresponding indices of each of its points
        tri_arr[mesh_tris[i][0]].append(i)
        tri_arr[mesh_tris[i][1]].append(i)
        tri_arr[mesh_tris[i][2]].append(i)

    # Sys of eq.s array
    fea_arr = np.zeros([len(mesh_points), len(mesh_points) + 1])
    mesh_obj_arr = np.empty(len(mesh_points), dtype=object)
    idx = 0
    for pt in mesh_points:
        adj = np.array(adj_arr[idx]) # Nodes adjacent to current node
        tri = np.empty(shape=[len(tri_arr[idx]), 3]) # Triangles containing current node as a vertex
        for i in range(len(tri_arr[idx])):
            # For each triangle
            # Rearrange so that the current node is in index 0

            triplet = mesh_tris[tri_arr[idx][i]]
            if triplet[1]==idx:
                triplet = np.roll(triplet, -1)
            elif triplet[2] ==idx:
                triplet = np.roll(triplet, -2)
            tri[i] = np.array(triplet)

        if len(adj) > np.shape(tri)[0]:
            # Boundary point - more adj edges than triangles
            bdry_bool = True
            i = 0
            while tri[i, 1] in tri[:, 2]:
                # Triangles are organized such that points are anticlockwise
                # Were arranged above so that each triplet has the current node at idx 0, then the remaining two ccw
                # Thus, to find boundary, cycle through values at idx 1 until you find one that doesn't show up in idx 2
                # This is the node directly ccw from the open region
                # i is the location of this node in tri[:][1]
                i+=1
            head = tri[i][1] # Node number of point directly ccw from boundary region
            # adj has no order in relation to tri
            # define ccw_start as idx in adj of head
            ccw_start = np.where(head==adj)

            rel_coords = np.array([mesh_points[i] - pt for i in adj]) # Coords of adj nodes relative to current node
            angle = np.array([np.mod(np.arctan2(rel_coords[i][1], rel_coords[i][0]), np.pi * 2) for i in range(np.shape(rel_coords)[0])])
            angle = np.mod(angle - angle[ccw_start], 2 * np.pi)

            inds = np.argsort(angle)
            adj = np.array([adj[i] for i in inds]) # Rearrange adj so nodes are in ccw order
            rel_coords = np.array([rel_coords[i] for i in inds]) # Rearrange so relative coords are in ccw order
            angle = np.array([angle[i] for i in inds]) # Rearrange so angles are in ccw order

            angle_area = np.empty(len(angle))
            angle_area[0] = angle[1]/2 # First angle - 1/2 angle from this node to second node (doesn't include boundary)
            # Middle angles - 1/2 angle from previous node + 1/2 angle to next node
            angle_area[1:len(angle)-1]=[(angle[i+1]-angle[i-1])/2 for i in np.arange(1, len(angle)-1)]
            angle_area[-1] = (angle[-1]-angle[-2])/2 # Third angle - 1/2 angle from prev node to this node (w/o boundary)
            angle_area = angle_area / (2*np.pi) # Normalize angles, since these are used as "areas" across which heat transfers
            angle_area = np.append(angle_area, 1 - np.sum(angle_area))


        else:
            # Interior point - same number of adj edges and triangles
            bdry_bool = False
            head = tri[0][2]
            # adj has no order in relation to tri
            # No boundary to start from
            # Simply define "head" as whatever the first adj node in tri happens to be
            ccw_start = np.where(head==adj)

            rel_coords = np.array([mesh_points[i] - pt for i in adj]) # Coords of adj nodes relative to current node
            # Compute angle from x-axis [0, 2pi]
            angle = np.array([np.mod(np.arctan2(rel_coords[i][1], rel_coords[i][0]), np.pi * 2) for i in range(np.shape(rel_coords)[0])])
            # Compute angle from start node [0, 2pi]
            angle = np.mod(angle - angle[ccw_start], 2 * np.pi)

            # Sort nodes, coords, and angles in increasing order according to angle
            # Puts the nodes in order ccw from start
            inds = np.argsort(angle)
            adj = np.array([adj[i] for i in inds]) # Rearrange adj so nodes are in ccw order
            rel_coords = np.array([rel_coords[i] for i in inds]) # Rearrange so relative coords are in ccw order
            angle = np.array([angle[i] for i in inds]) # Rearrange so angles are in ccw order

            # Compute arc length for each adjacent node
            angle_area = np.empty(len(angle))
            angle_area[0] = (angle[1]+(2*np.pi - angle[-1]))/2 # First angle - 1/2 angle from this node to second node (doesn't include boundary)
            # Middle angles - 1/2 angle from previous node + 1/2 angle to next node
            angle_area[1:len(angle)-1]=[(angle[i+1]-angle[i-1])/2 for i in np.arange(1, len(angle)-1)]
            angle_area[-1] = (2*np.pi - angle[-2])/2 # Third angle - 1/2 angle from prev node to this node (w/o boundary)
            angle_area = angle_area / (2*np.pi) # Normalize angles, since these are used as "areas" across which heat transfers


        eq = np.zeros(len(mesh_points)) #Heat balance equation for this point
        for i in range(len(adj)):
            t = np.linalg.norm(rel_coords[i])
            eq[adj[i]] += k/t * angle_area[i]
            eq[idx] -= k/t * angle_area[i]

        fea_arr[idx, :-1] = eq


        point = MeshPoint(mesh, idx, pt, adj, rel_coords, angle_area, tri, eq, bdry_bool)
        mesh_obj_arr[idx] = point

        idx += 1
    '''
    fea_arr = np.zeros([len(mesh_points), len(mesh_points) + 1])
    for ptpair in mesh_faces:
        # Conductive heat transfer between each pair of points
        p1 = ptpair[0] # Separate pair of points
        p2 = ptpair[1]

        # Compute Euclidean distance between points
        t = np.linalg.norm(mesh_points[p1] - mesh_points[p2])

        fea_arr[p1, p1] -= k_OW / t
        fea_arr[p1, p2] += k_OW / t
        fea_arr[p2, p1] += k_OW / t
        fea_arr[p2, p2] -= k_OW / t
    '''
    # Add effect of boundary conditions
    '''
    Boundary condition key:
    bc = [type, coords, value, dtype=object]
    type:
        0 - set temperature
        1 - set convection
    coords:
        [start_point, end_point]
        where start_point and end_point must be tuples that exist in mesh_points
    value
        temperature: input temperature of points in kelvin
        convection: [conv_coeff, amb_temp]
    '''

    for bc in np.arange(np.shape(boundary_conds)[1]-1):
        coords = np.array(boundary_conds[bc][1])
        coords = np.where((mesh_points == coords[:, None]).all(-1))[1]

        if boundary_conds[bc][0]==1: # Convection

            # Get coordinates of points to apply boundary condition, account for looping back to start if needed
            if coords[0] > coords[1]:
                points = np.arange(coords[0], len(bdry))
                points = np.append(points, 0)
            else:
                points = np.arange(coords[0],coords[1]+1)

            for i in np.arange(np.shape(points)[0]):
                point = mesh_obj_arr[points[i]]
                fea_arr[points[i], points[i]] -= boundary_conds[bc][2][0] * point.areas[-1]
                fea_arr[points[i], -1] -= boundary_conds[bc][2][0]*boundary_conds[bc][2][1] * point.areas[-1]

        # Keep this at end so that it can override other boundary conditions
        if boundary_conds[bc][0] == 0:
            if coords[0] > coords[1]:
                points = np.arange(coords[0], len(bdry))
                points = np.append(points, 0)
            else:
                points = np.arange(coords[0],coords[1]+1)

            for i in np.arange(np.shape(points)[0]):
                coeffs = fea_arr[:, points[i]]
                fea_arr[:, -1] -= coeffs * boundary_conds[bc][2]
                fea_arr[points[i],:] = np.zeros(np.shape(fea_arr)[1])
                fea_arr[:, points[i]] = np.zeros(np.shape(fea_arr)[0])
                fea_arr[points[i], points[i]] = 1
                fea_arr[points[i], -1] = boundary_conds[bc][2]

    print(np.shape(fea_arr))

    LU, piv = la.lu_factor(fea_arr[:, :-1])
    temps = la.lu_solve((LU,piv), fea_arr[:, -1])


    x = []
    y = []
    for i in np.arange(np.shape(mesh_points)[0]):
        x = np.append(x, mesh_points[i, 0])
        y = np.append(y, mesh_points[i, 1])
    '''
    x = x.astype('float')
    y = y.astype('float')
    temps = temps.astype('float')
    plt.figure()
    plt.tricontourf(x,y,np.asarray(mesh_tris),temps, levels=10, cmap='jet')
    plt.colorbar()
    '''
    return (mesh, mesh_points, temps)

'''
res = .00005 # Resolution aka max distance between edge points
corner_points = [(0, 0), (0, .003), (.001, .003), (.001, .001), (.002, .001), (.002, 0)]
boundary_conds = np.array([
    [0, [corner_points[5], corner_points[0]], 350],
    [1, [corner_points[2], corner_points[4]], [1000, 300]],
], dtype=object)
k_OW = 250
fea_solver(corner_points, res, boundary_conds, k_OW)
plt.show()
'''