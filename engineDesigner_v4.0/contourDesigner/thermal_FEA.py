import matplotlib.pyplot as plt
import meshpy.triangle as triangle
import numpy as np
import scipy.linalg as la


def generate_mesh(corner_points, res):
    def round_trip_connect(start, end):
        return [(i, i + 1) for i in range(start, end)] + [(end, start)]

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
    mesh_info = triangle.MeshInfo() # Generate info for mesh
    mesh_info.set_points(bdry)  # Mesh boundary points
    mesh_info.set_facets(bdry_faces) # Mesh boundary edges

    mesh = triangle.build(mesh_info, generate_faces=True, volume_constraints=True,
                          max_volume=res ** 2) # Build mesh - note max triangle area set to res^2 (e.g. a square of side length res)

    mesh_points = np.array(mesh.points)
    mesh_faces = np.array(mesh.faces)
    mesh_tris = np.array(mesh.elements)
    '''
    plt.figure()
    plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris)
    plt.xlabel("x")
    plt.ylabel("y")
    '''

    return (mesh, mesh_points, mesh_faces, mesh_tris, bdry, bdry_faces)

def fea_solver(corner_points, res, boundary_conds, k):

    #Generate mesh from boundary
    (mesh, mesh_points, mesh_faces, mesh_tris, bdry, bdry_facets) = generate_mesh(corner_points, res)

    # Sys of eq.s array
    fea_arr = np.zeros([len(mesh_points), len(mesh_points) + 1])
    for ptpair in mesh_faces:
        # Conductive heat transfer between each pair of points
        p1 = ptpair[0] # Separate pair of points
        p2 = ptpair[1]

        # Compute Euclidean distance between points
        t = np.linalg.norm(mesh_points[p1] - mesh_points[p2])

        fea_arr[p1, p1] -= k / t
        fea_arr[p1, p2] += k / t
        fea_arr[p2, p1] += k / t
        fea_arr[p2, p2] -= k / t

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

        if boundary_conds[bc][0]==1:

            if coords[0] > coords[1]:
                points = np.arange(coords[0], len(bdry))
                points = np.append(points, 0)
            else:
                points = np.arange(coords[0],coords[1]+1)

            for i in np.arange(np.shape(points)[0]):
                fea_arr[points[i], points[i]] -= boundary_conds[bc][2][0]*3/2
                fea_arr[points[i], -1] -= boundary_conds[bc][2][0]*boundary_conds[bc][2][1]*3/2

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

    '''
    fea_arr = sympy.Matrix(fea_arr)
    fea_arr_sol = fea_arr.rref()[0]
    fea_arr_sol = np.array(fea_arr_sol)
    print(fea_arr_sol)
    '''
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
    plt.tricontourf(x,y,np.asarray(mesh_tris),temps)
    '''
    return (mesh, mesh_points, temps)

'''
res = .00005  # Resolution aka max distance between edge points
corner_points = [(0, 0), (0, .003), (.001, .003), (.001, .001), (.002, .001), (.002, 0)]
boundary_conds = np.array([
    [0, [corner_points[5], corner_points[0]], 350],
    [1, [corner_points[2], corner_points[4]], [1000, 300]],
], dtype=object)
k_OW = 50
fea_solver(corner_points, res, boundary_conds, k_OW)
plt.show()
'''