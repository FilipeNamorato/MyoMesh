import dolfin as df
import ldrb
import meshio


def solve_laplace(mesh, boundary_markers, boundary_values, ldrb_markers):
    V = df.FunctionSpace(mesh, 'P', 1)

    u_rv, u_lv, u_epi = boundary_values

    bc1 = df.DirichletBC(V, u_rv, boundary_markers, ldrb_markers["rv"]) 
    bc2 = df.DirichletBC(V, u_lv, boundary_markers, ldrb_markers["lv"])
    bc3 = df.DirichletBC(V, u_epi, boundary_markers, ldrb_markers["epi"])

    bcs=[bc1, bc2 ,bc3]

    ds = df.Measure('ds', domain=mesh, subdomain_data=boundary_markers)
    dx = df.Measure('dx', domain=mesh)

    # Define variational problem
    u = df.TrialFunction(V)
    v = df.TestFunction(V)
    f = df.Constant(0.0)   
    a = df.dot(df.grad(u), df.grad(v))*dx  
    L = f*v*dx

    # Compute solution
    u = df.Function(V)
    df.solve(a == L, u, bcs, solver_parameters=dict(linear_solver='gmres', preconditioner='hypre_amg')) 

    return u


def request_functions(meshname, aux_alpha_endo_lv, aux_alpha_epi_lv, aux_beta_endo_lv, 
                    aux_beta_epi_lv, aux_alpha_endo_sept, aux_alpha_epi_sept, 
                    aux_beta_endo_sept, aux_beta_epi_sept, aux_alpha_endo_rv, 
                    aux_alpha_epi_rv, aux_beta_endo_rv, aux_beta_epi_rv):

    #mesh reading converted by dolfin-convert
    mesh = df.Mesh(meshname + '.xml')
    materials = df.MeshFunction("size_t", mesh, meshname + '_physical_region.xml')
    ffun = df.MeshFunction("size_t", mesh, meshname + '_facet_region.xml')

    V0 = df.FunctionSpace(mesh, "DG", 0) 
    tecido = df.Function(V0) 
    tecido.vector()[:] = materials.array()==2
    tecido.rename("tecido", "tecido")

    ldrb_markers = {
        "base": 10,
        "lv": 30,
        "epi": 40,
        "rv": 20
    }

    # Choose space for the fiber fields
    # This is a string on the form {family}_{degree}
    fiber_space = "DG_0"

    #--------
    # Create field to define the action potential fenotype
    # Solve Laplace problems with different boundary conditions
    # u=1 on epicardium
    phi_epi = solve_laplace(mesh, ffun, [0, 0, 1], ldrb_markers)
    # u=1 on LV endocardium
    phi_lv = solve_laplace(mesh, ffun, [0, 1, 0], ldrb_markers)
    # u=1 on RV endocardium
    phi_rv = solve_laplace(mesh, ffun, [1, 0, 0], ldrb_markers)

    # Compute field with Laplace solutions
    V = df.FunctionSpace(mesh, 'Lagrange', 1)
    u = df.Function(V)
    u.interpolate(df.Expression('-(epi + 2*rv*lv/(rv+lv) ) + 1', epi=phi_epi, rv=phi_rv, lv=phi_lv, degree=1))

    bc3 = df.DirichletBC(V, 0, ffun, ldrb_markers["epi"])
    bc3.apply(u.vector())

    u.rename("fenotipo","fenotipo")
    #--------

    # Compute the microstructure
    fiber, sheet, sheet_normal = ldrb.dolfin_ldrb(
        mesh=mesh,
        fiber_space=fiber_space,
        ffun=ffun,
        markers=ldrb_markers,
        alpha_endo_lv=aux_alpha_endo_lv,  # Fiber angle on the LV endocardium
        alpha_epi_lv=aux_alpha_epi_lv,  # Fiber angle on the LV epicardium
        beta_endo_lv=aux_beta_endo_lv,  # Sheet angle on the LV endocardium
        beta_epi_lv=aux_beta_epi_lv,  # Sheet angle on the LV epicardium
        alpha_endo_sept=aux_alpha_endo_sept,  # Fiber angle on the Septum endocardium
        alpha_epi_sept=aux_alpha_epi_sept,  # Fiber angle on the Septum epicardium
        beta_endo_sept=aux_beta_endo_sept,  # Sheet angle on the Septum endocardium
        beta_epi_sept=aux_beta_epi_sept,  # Sheet angle on the Septum epicardium
        alpha_endo_rv=aux_alpha_endo_rv,  # Fiber angle on the RV endocardium
        alpha_epi_rv=aux_alpha_epi_rv,  # Fiber angle on the RV epicardium
        beta_endo_rv=aux_beta_endo_rv,  # Sheet angle on the RV endocardium
        beta_epi_rv=aux_beta_epi_rv,
    )

    fiber.rename("f_0","f_0")
    sheet.rename("s_0","s_0")
    sheet_normal.rename("n_0","n_0")

    print("Saving...")

    with df.XDMFFile(mesh.mpi_comm(), meshname + ".xdmf") as xdmf:
        xdmf.parameters.update(
        {
            "functions_share_mesh": True,
            "rewrite_function_mesh": False
        })
        xdmf.write(mesh)
        xdmf.write(fiber, 0)
        xdmf.write(sheet, 0)
        xdmf.write(sheet_normal,0)
        xdmf.write(tecido, 0)
        xdmf.write(u, 0)


    convert_xdmf_to_vtu(meshname)

    print("Done.")

def convert_xdmf_to_vtu(meshname):

    print("Converting .xdmf to .vtu")
    filename = meshname+".xdmf"
    t, point_data, cell_data = None, None, None

    with meshio.xdmf.TimeSeriesReader(filename) as reader:
        points, cells = reader.read_points_cells()
        t, point_data, cell_data = reader.read_data(0)

    mesh = meshio.Mesh(points, cells, point_data=point_data, cell_data=cell_data,)
    mesh.write(meshname+".vtu", file_format="vtu",  )