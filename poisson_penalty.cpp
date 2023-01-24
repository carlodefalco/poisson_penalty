#include <bim_sparse.h>
#include <limits>
#include <mumps_class.h>
#include <quad_operators.h>
#include <tmesh.h>
#include <vector>


#include "poisson_penalty.h"



static double
NACA (double x) {
  return  5 * h2 *
    (0.2969 * std::sqrt(x) + x * ( -.1260 + x * (- .3516 + x * (.2843 - .1015 * x))));
}

static double
distfun_airfoil (double x, double y){
  if ((x-d)/a <= 0) {
    return - y*y - std::pow((x-d)/a, 2);
  } else if ((x-d)/a >= 1) {
    return - y*y - std::pow((x-d)/a, 2);
  } else {
    return std::pow(NACA((x-d)/a), 2) - y*y;
  }
  return 0.;
}

static double
distfun_source (double x, double y){
  return std::pow (1 - (std::pow((x-(-d-lc/2))/lc*2,2) + std::pow(y/h1,2)), 5.) ;
}



static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return 1; }

static int
contact_refinement (tmesh::quadrant_iterator q)
{
  double x, y;
  for (int ii = 0; ii < 4; ++ii) {
    x = q->p(0, ii);
    y = q->p(1, ii);
    
    if (y<=h1 && x>=-d-lc && x <=-d) {
      return 1;
    }

    if (y<=h2 && x>=d && x <=d+a) {
      return 1;
    }
  }
  
  return 0;
}

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh                 tmsh;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  /*
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < 3; ++cycle)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }
  */
  
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < 3; ++cycle)
    {
      tmsh.set_refine_marker (contact_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }

  tmsh.vtk_export ("poisson_penalty");
  
  std::vector<tmesh::idx_t> nnodes;
  std::vector<double> h_step;
  // std::vector<double> error;
  
  for (int adapt = 0; adapt < refine_steps; ++adapt)
    {
      std::cout << "*** Step " << adapt << " ***" << std::endl;
      
      // Compute coefficients.      
      std::vector<double> alpha(tmsh.num_global_nodes (), kG);
      std::vector<double> beta(tmsh.num_global_nodes (), 0.0);
      std::vector<double> psi(tmsh.num_global_nodes (), 0.0);
      
      std::vector<double> f(tmsh.num_local_quadrants (), 1);
      std::vector<double> g(tmsh.num_global_nodes (), 0.0);
      std::vector<double> df(tmsh.num_global_nodes (), 0.0);
      
      double x = 0, y = 0;
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          for (int ii = 0; ii < 4; ++ii)
            {
              x = quadrant->p(0, ii);
              y = quadrant->p(1, ii);
                
              if (! quadrant->is_hanging (ii) ) {
		if (distfun_airfoil (x, y) >= 0.) {
		  alpha[quadrant->gt(ii)] = kS;
		  beta[quadrant->gt(ii)] = 1.;
		  g[quadrant->gt(ii)] = 10.e3;
		}
		if (distfun_source (x, y) >= 0.) {
		  alpha[quadrant->gt(ii)] = kS;
		  beta[quadrant->gt(ii)] = 1.;
		}
		df[quadrant->gt(ii)] = distfun_airfoil (x, y) * distfun_source (x, y);
	      }
            }
        }
      
      // Assemble system matrix and right-hand side.
      sparse_matrix A;
      A.resize(tmsh.num_global_nodes());
      
      // Reduce coefficients.
      std::vector<double> global_alpha(tmsh.num_global_nodes(), 0);
      MPI_Allreduce(alpha.data(), global_alpha.data(), alpha.size(),
                    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      
      bim2a_advection_eafe_diffusion (tmsh, global_alpha, psi, A);
      bim2a_reaction (tmsh, f, beta, A);
      
      std::vector<double> rhs(tmsh.num_global_nodes (), 0);
      bim2a_rhs (tmsh, f, g, rhs);
      
      // Set boundary conditions.

      /*
      func u_ex =
        [] (double x, double y)
        {
	  return 0.;
        };
      */
      
      // dirichlet_bcs bcs;
      // for (int i = 0; i < 4; ++i)
      //   bcs.push_back (std::make_tuple(0, i, u_ex));
      
      // bim2a_dirichlet_bc (tmsh, bcs, A, rhs);
      
      // Solve problem.
      std::cout << "Solving linear system.";
      
      mumps mumps_solver;
      
      std::vector<double> vals;
      std::vector<int> irow, jcol;
      
      A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
      
      mumps_solver.set_lhs_distributed ();
      mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
      mumps_solver.set_distributed_lhs_data (vals);
      
      // Reduce rhs (so that rank 0 has the actual rhs).
      std::vector<double> global_rhs(tmsh.num_global_nodes(), 0);
      MPI_Reduce(rhs.data(), global_rhs.data(), rhs.size(),
                 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      
      if (rank == 0)
        mumps_solver.set_rhs (global_rhs);
      
      // Solve.
      mumps_solver.analyze ();
      mumps_solver.factorize ();
      mumps_solver.solve ();
      mumps_solver.cleanup ();
      
      // Export solution.
      MPI_Bcast(global_rhs.data(), global_rhs.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

      tmsh.octbin_export ((std::string("poisson_penalty_u_")
                           + std::to_string(adapt)).c_str(), global_rhs);

      tmsh.octbin_export ((std::string("poisson_penalty_beta_")
                           + std::to_string(adapt)).c_str(), beta);

      tmsh.octbin_export ((std::string("poisson_penalty_g_")
                           + std::to_string(adapt)).c_str(), g);

      tmsh.octbin_export ((std::string("poisson_penalty_distfun_")
			   + std::to_string(adapt)).c_str(), df);
      
      //std::vector<double> uex(tmsh.num_global_nodes(), 0);
      
      // for (auto quadrant = tmsh.begin_quadrant_sweep ();
      //      quadrant != tmsh.end_quadrant_sweep ();
      //      ++quadrant)
      //   for (int i = 0; i < 4; ++i)
      //     uex[quadrant->gt(i)] = u_ex(quadrant->p(0, i), quadrant->p(1, i));
      
      // tmsh.octbin_export ((std::string("poisson_penalty_uex_")
      //                      + std::to_string(adapt)).c_str(), uex);
      
      std::cout << " Done." << std::endl;
      
      // Compute reconstructed gradient.
      std::cout << "Computing reconstructed gradient and estimator.";
      
      active_fun regionG = [] (tmesh::quadrant_iterator q)
        { return (distfun_airfoil (q->centroid(0), q->centroid(1)) < 0.)
	  && (distfun_source (q->centroid(0), q->centroid(1)) < 0.); };
      
      active_fun regionS = [] (tmesh::quadrant_iterator q)
        { return (distfun_airfoil (q->centroid(0), q->centroid(1)) >= 0.)
	  || (distfun_source (q->centroid(0), q->centroid(1)) >= 0.); };
      
      gradient<std::vector<double>> du0 = bim2c_quadtree_pde_recovered_gradient(tmsh, global_rhs, regionG);
      gradient<std::vector<double>> du1 = bim2c_quadtree_pde_recovered_gradient(tmsh, global_rhs, regionS);
      
      q2_vec u_star0 = bim2c_quadtree_pde_recovered_solution(tmsh, global_rhs, du0);
      q2_vec u_star1 = bim2c_quadtree_pde_recovered_solution(tmsh, global_rhs, du1);
      
      tmsh.octbin_export ((std::string("poisson_penalty_du0_x_")
                           + std::to_string(adapt)).c_str(), du0.first);
      tmsh.octbin_export ((std::string("poisson_penalty_du0_y_")
                           + std::to_string(adapt)).c_str(), du0.second);
      
      tmsh.octbin_export ((std::string("poisson_penalty_du1_x_")
                           + std::to_string(adapt)).c_str(), du1.first);
      tmsh.octbin_export ((std::string("poisson_penalty_du1_y_")
                           + std::to_string(adapt)).c_str(), du1.second);
      
      auto estimator = [& du0, & du1, & global_rhs] (tmesh::quadrant_iterator q)
        {
          if ((distfun_airfoil (q->centroid(0), q->centroid(1)) < 0.)
	      && (distfun_source (q->centroid(0), q->centroid(1)) < 0.))
            return estimator_grad (q, du0, global_rhs);
          else
            return estimator_grad (q, du1, global_rhs);
        };
      
      std::cout << " Done." << std::endl;
      
      // Compute h and error.
      double hx = 0, hy = 0,
        h = std::numeric_limits<double>::max (),
        global_h = 0;
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          hx = quadrant->p(0, 1) - quadrant->p(0, 0);
          hy = quadrant->p(1, 2) - quadrant->p(1, 0);
          
          h = std::min(h, std::sqrt(hx*hx + hy*hy));
          
        }
      
      MPI_Reduce(&h, &global_h, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);
      
      nnodes.push_back (tmsh.num_global_nodes ());
      h_step.push_back (global_h);
      
      std::cout << " Done." << std::endl;
      
      if (tmsh.num_global_nodes () >= 1e6)
        break;
      
      // Refine.
      tmsh.set_metrics_marker (estimator, 1e-4, 7, 1, 1);
      tmsh.metrics_refine (1e4);
      
      tmsh.vtk_export ((std::string("poisson_penalty_newmesh_")
                        + std::to_string(adapt)).c_str());
    }
  
  if (rank == 0)
    for (unsigned step = 0; step < nnodes.size(); ++step)
      std::cout << "Step " << step << ", #nodes: "
                << nnodes[step] << ", h: "
                << h_step[step] 
		<< std::endl;
  
  MPI_Finalize ();
  
  return 0;
}
