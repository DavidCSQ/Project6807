#pragma once

#define IGL_VIEWER_VIEWER_QUIET 1
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/ViewerPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/lbs_matrix.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject_in_mesh.h>
#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/decimate.h>

#include "handles.hpp"
#include "mesh_util.hpp"
#include "nearest_neighbor_weights.hpp"
#include "bounded_biharmonic_weights.hpp"
#include "linear_weights.hpp"

class MakeItStandBasicUI : public igl::opengl::glfw::ViewerPlugin {
public:
	MakeItStandBasicUI() :
		create_handle_mode(false),
		weight_type(0),
		menu()
	{
		// Nothing to do
	}

	void init(igl::opengl::glfw::Viewer* _viewer) {
		igl::opengl::glfw::ViewerPlugin::init(_viewer);
		viewer->core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL);

		viewer->plugins.push_back(&menu);
		menu.callback_draw_viewer_menu = [&]() {draw_menu(); };

	}

	void draw_menu() {

		// Helper for making checkboxes
		auto make_checkbox = [&](const char* label, unsigned int& option)
		{
			return ImGui::Checkbox(label,
				[&]() { return viewer->core().is_set(option); },
				[&](bool value) { return viewer->core().set(option, value); }
			);
		};

		if (ImGui::Button("Load Triangle Mesh")) {
			auto filename = igl::file_dialog_open();
			if (filename.size() > 0) {
				Eigen::MatrixXd V_temp, V_temp_unique;
				Eigen::MatrixXi F_temp;
				Eigen::VectorXi SVI, SVJ;
				igl::read_triangle_mesh(filename, V_temp, F_temp);
				igl::remove_duplicate_vertices(V_temp, 0, V_temp_unique, SVI, SVJ);
				std::for_each(F_temp.data(), F_temp.data() + F_temp.size(), [&SVJ](int& f) {f = SVJ(f); });
				V_temp = V_temp_unique;
				std::cout << "Triangle Mesh has faces: " << F_temp.rows() << std::endl;

				if (F_temp.rows() > MAX_FACES) {
					std::cout << "Mesh too large, decimating." << std::endl;
					Eigen::MatrixXd V_dec_temp;
					Eigen::MatrixXi F_dec_temp;
					Eigen::VectorXi J, I;

					std::cout << "Size before = " << V_temp.rows() << " x " << V_temp.cols() << std::endl;

					if (igl::decimate(V_temp, F_temp, MAX_FACES, V_dec_temp, F_dec_temp, J, I)) {
						std::cout << "Size after = " << V_dec_temp.rows() << " x " << V_dec_temp.cols() << std::endl;
						V_temp = V_dec_temp;
						F_temp = F_dec_temp;
					}
					else {
						std::cout << "Decimation failed" << std::endl;
					}
				}

				std::cout << "Normal vertices: " << V_temp.rows() << std::endl;
				if (igl::copyleft::tetgen::tetrahedralize(V_temp, F_temp, "pq1.414", V, T, F) == 0) {
					F.col(0).swap(F.col(2));

					viewer->data().clear();
					viewer->data().set_mesh(V, F);
					viewer->core().align_camera_center(viewer->data().V, viewer->data().F);
					handles = WeightHandles();
					create_handle_mode = 0;
				}
				else {
					std::cout << "Failed to tetrahedralize mesh." << std::endl;
				}
				std::cerr << "number of vertices = " << V.rows() << std::endl;
				std::cerr << "number of faces = " << T.rows() << std::endl;

				set_support_vertices();
			}
		}

		// Only Allow Mesh Manipulations if we have a mesh
		if (V.rows() > 0) {

			ImGui::NewLine();

			if (ImGui::Button("Load Handles")) {
				auto filename = igl::file_dialog_open();
				if (filename.size() > 0) {
					handles.load_handle_file(filename);
					draw_handles();
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("Save Handles")) {
				auto filename = igl::file_dialog_save();
				if (filename.size() > 0) {
					handles.save_handle_file(filename);
				}
			}

			ImGui::NewLine();
			ImGui::Text("Handle Mode");
			if (ImGui::RadioButton("Create", &create_handle_mode, 0)) {
				handles.reset_handle_transformations();
				viewer->data().set_vertices(V);
				draw_handles();
				//fixed = false;
			}
			// Only allow handle manipulation when there are handles to manipulate
			if (handles.positions().rows() > 0) {
				ImGui::SameLine();
				if (ImGui::RadioButton("Manipulate", &create_handle_mode, 1)) {
					started_mesh_manipulation = true;
					if (create_handle_mode == 1) {
						std::cout << "Recomputing Weights - UI may hang during this." << std::endl;
						compute_weights();
					}
				}
			}

			ImGui::NewLine();
			
			if (!started_mesh_manipulation) {
				ImGui::Text("Support Point Threshold:");
				ImGui::DragFloat("Threshold", &support_threshold, 0.1f, 0.001 * (max_vertex_height - min_vertex_height), 0.1 * (max_vertex_height - min_vertex_height));
				if (ImGui::Button("Recompute Support")) {
					set_support_vertices();
					draw_handles();
					draw_support_polygon();
				}

				ImGui::NewLine();
			}

			if (ImGui::Button("Compute Center Of Gravity")) {
				compute_cog();
				draw_handles();
				draw_support_polygon();
			}

			if (cog_computed) {
				ImGui::Text("Will it stand?");
				std::string will_stand;
				if (will_it_stand(c)) will_stand = "Yes";
				else will_stand = "No";
				ImGui::Text(will_stand.c_str());
			}
			if (ImGui::Button("Optimize current view")) {
				optimize();
			}

			if (ImGui::Button("Fix all the created handles")) {
				n_handles = handles.positions().rows();
			}
			ImGui::NewLine();
			if (ImGui::CollapsingHeader("Viewer Options")) {
				make_checkbox("Wireframe", viewer->data().show_lines);
				make_checkbox("Fill", viewer->data().show_faces);
			}
		}
	}
    void optimize(){//assume in manipulate mode
        int maxiter = 10;
        for(int i = 0;i<maxiter;i++){
		    //compute_weights();
            compute_w_opt();
            Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Upper> solver;
            solver.setMaxIterations(1000000); //just a large number
            Eigen::MatrixXd K = dch.block(0,n_handles*3,3,(handles.positions().rows()-n_handles)*3);
            Eigen::Vector3d fe = Eigen::Vector3d(sc[0]-c[0],sc[1]-c[1],0);
			if (fe.norm()<1e-6) break;
			Eigen::VectorXd u = solver.compute(K).solve(fe);
			Eigen::MatrixXd H;
			handles.transfored_handles(H);
			std::cout<<"ftest\n" <<K*u+c<<"\n current pos\n"<<H.row(n_handles)<<"\n u\n" << u << std::endl;
			//Eigen::VectorXd u = Eigen::Vector3d(0.03,0.03,0.03);
			//Eigen::Vector3d ftest = K*u + c;
            //std::cout << "Time to solve" << "\n K\n" << K << "\n fe\n" << fe << "\n c\n" << c <<"\n u\n" << u << std::endl;
		    //std::cout << "ftest\n" << ftest << std::endl;
			for(int j = n_handles;j<handles.positions().rows();j++){
				Eigen::RowVector3d position = handles.transform().block(4 * j + 3, 0, 1, 3);
			    Eigen::Vector3d pos = Eigen::Vector3d(position[0]+handles.positions().row(j)[0]+u[(j-n_handles)*3],position[1]+handles.positions().row(j)[1]+u[(j-n_handles)*3+1],position[2]+handles.positions().row(j)[2]+u[(j-n_handles)*3+2]);
                handles.move_handle(
				    j,
				    pos);
			    std::cout<<"newpos"<<pos;
			}
			//std::cout<<"nhandles"<<n_handles<<"max"<<handles.positions().rows();
			viewer->data().set_vertices(lbs_mat * handles.transform());
			draw_handles();
			compute_cog();
			//if (fe.norm()<1e-6) break;
			
		}
		double ddm = 0;
		Eigen::Vector3d ddc;
		ddc.setZero();
		for(int i = 0; i < viewer->data().V.rows();i++){
            ddm += (dm[i*3]*(viewer->data().V.row(i)[0]-lbs_mat(i,0))+dm[i*3+1]*(viewer->data().V.row(i)[1]-lbs_mat(i,1))+dm[i*3+2]*(viewer->data().V.row(i)[2]-lbs_mat(i,2)))/18;
			ddc += (dcf.col(i*3)*(viewer->data().V.row(i)[0]-lbs_mat(i,0))+dcf.col(i*3+1)*(viewer->data().V.row(i)[1]-lbs_mat(i,1))+dcf.col(i*3+2)*(viewer->data().V.row(i)[2]-lbs_mat(i,2)));
		}
		//std::cout<< "ddm  " << ddm <<"  ddc  "<<ddc<<std::endl;
	}
	Eigen::Vector3d direction_to_support_polyggon() {
		// Need to compute convex hull of a base plane
		double bottom = V.col(2).minCoeff();

		for (int f = 0; f < F.rows(); f++) {
			auto vinds = F.row(f);

			const Eigen::Vector3d vi = deformed_V.row(vinds[0]);
			const Eigen::Vector3d vj = deformed_V.row(vinds[1]);
			const Eigen::Vector3d vk = deformed_V.row(vinds[2]);
			
			// We assume the user chose the support polygon to be the base plane of the mesh
		}
	}
    void compute_w_opt(){
		compute_cog();
		dch.resize(3,3*handles.positions().rows());
		dch.setZero();
        for(int f = 0;f<handles.positions().rows();f++){
			for(int l = 0;l<lbs_mat.rows();l++){
                dch(0,f*3)+=dcf(0,3*l)*lbs_mat(l,f*4+3);
				dch(0,f*3+1)+=dcf(0,3*l+1)*lbs_mat(l,f*4+3);
				dch(0,f*3+2)+=dcf(0,3*l+2)*lbs_mat(l,f*4+3);
				dch(1,f*3+1)+=dcf(1,3*l+1)*lbs_mat(l,f*4+3);
				dch(1,f*3)+=dcf(1,3*l)*lbs_mat(l,f*4+3);
				dch(1,f*3+2)+=dcf(1,3*l+2)*lbs_mat(l,f*4+3);
				dch(2,f*3+2)+=dcf(2,3*l+2)*lbs_mat(l,f*4+3);
				dch(2,f*3)+=dcf(2,3*l)*lbs_mat(l,f*4+3);
				dch(2,f*3+1)+=dcf(2,3*l+1)*lbs_mat(l,f*4+3);
			}
		}
	}
	void compute_weights() {
		// Only compute weights if there are weights to compute!
		Eigen::MatrixXd W;
		Eigen::VectorXi b;
		Eigen::MatrixXd bc;
		compute_cog();
		handles.boundary_conditions(V, T, b, bc);
		bounded_biharmonic_weights(V, T, b, bc, W);
		igl::lbs_matrix(V, W, lbs_mat);
	}

	Eigen::Vector3d cross_prod(Eigen::Vector3d a, Eigen::Vector3d b) {
		double element1 = a[1] * b[2] - a[2] * b[1];
		double element2 = a[2] * b[0] - a[0] * b[2];
		double element3 = a[0] * b[1] - a[1] * b[0];
		return { element1, element2, element3 };
	}

	void compute_cog() {
		// The authors made a lot of mistakes in their code, but I am sticking to their approach for now
		// Eigen::Matrix<double, 1, Eigen::Dynamic> mList;
		// Eigen::Matrix<double, 3, Eigen::Dynamic> cList;
		Eigen::Matrix<double, 3, Eigen::Dynamic> dmList;
		Eigen::Matrix<double, 9, Eigen::Dynamic> dcList;

		double m = 0.f;
		c.setZero();

		dmList.resize(3, 3 * F.rows());

		// We will compute derivatives at the same time we do mass and cog
		dm.resize(1, 3 * V.rows());
		dc.resize(3, 3 * V.rows());
		dm.setZero();
		dc.setZero();
		dcf.resize(3, 3 * V.rows());
		dcf.setZero();

		deformed_V = V;
		if (handles.positions().rows() > 0) {
			deformed_V = viewer->data().V;
		}

		std::cout << "V sample: " << deformed_V.row(100) << std::endl;

		// Iterate over faces to calculate contributions to mass and center of gravity (and derivatives) using divergence theorem
		for (int f = 0; f < viewer->data().F.rows(); f++) {
			// Get three vertices of face
			auto vinds = viewer->data().F.row(f);

			const Eigen::Vector3d vi = viewer->data().V.row(vinds[0]);
			const Eigen::Vector3d vj = viewer->data().V.row(vinds[1]);
			const Eigen::Vector3d vk = viewer->data().V.row(vinds[2]);

			
			Eigen::Vector3d e1 = vj - vi, e2 = vi - vk, e3 = vk - vj;
			// auto normal = -e1.cross(e2); This is broken
			auto normal = cross_prod(-e1, e2);

			// std::cout << "Face " << f << std::endl;
			// std::cout << "Vertex 0: " << std::endl << vi << std::endl;
			// std::cout << "Vertex 1: " << std::endl << vj << std::endl;
			// std::cout << "Vertex 2: " << std::endl << vk << std::endl;
			// std::cout << "Normal: " << std::endl << normal << std::endl;

			Eigen::Vector3d vsum = vi + vj + vk;
			Eigen::Vector3d g = vsum.cwiseProduct(vsum) - (
				vi.cwiseProduct(vj) + vj.cwiseProduct(vk) + vk.cwiseProduct(vi));
			// Eigen::Vector3d g = vi.cwiseProduct(vi) + vi.cwiseProduct(vj) + 
			// 	vj.cwiseProduct(vj) + vj.cwiseProduct(vk) + vk.cwiseProduct(vk) + vk.cwiseProduct(vi);

			// mass
			m += (vsum.dot(normal)); // This is the correct way according to divergence theorem
			// m += (vsum[0] * normal[0]); // This is the incorrect way the authors do it

			// center of mass
			c += g.cwiseProduct(normal);

			// derivative of mass contributions
			dm[3 * vinds[0]] += normal[0]-vk[1]*vsum[2]+vk[2]*vsum[1]-vj[2]*vsum[1]+vj[1]*vsum[2];
			dm[3 * vinds[0] + 1] += normal[1]-vk[2]*vsum[0]+vk[0]*vsum[2]-vj[0]*vsum[2]+vj[2]*vsum[0];
			dm[3 * vinds[0] + 2] += normal[2]-vk[0]*vsum[1]+vk[1]*vsum[0]-vj[1]*vsum[0]+vj[0]*vsum[1];

			dm[3 * vinds[1]] += normal[0]-vi[1]*vsum[2]+vi[2]*vsum[1]-vk[2]*vsum[1]+vk[1]*vsum[2];
			dm[3 * vinds[1] + 1] += normal[1]-vi[2]*vsum[0]+vi[0]*vsum[2]-vk[0]*vsum[2]+vk[2]*vsum[0];
			dm[3 * vinds[1] + 2] += normal[2]-vi[0]*vsum[1]+vi[1]*vsum[0]-vk[1]*vsum[0]+vk[0]*vsum[1];

			dm[3 * vinds[2]] += normal[0]-vj[1]*vsum[2]+vj[2]*vsum[1]-vi[2]*vsum[1]+vi[1]*vsum[2];
			dm[3 * vinds[2] + 1] += normal[1]-vj[2]*vsum[0]+vj[0]*vsum[2]-vi[0]*vsum[2]+vi[2]*vsum[0];
			dm[3 * vinds[2] + 2] += normal[2]-vj[0]*vsum[1]+vj[1]*vsum[0]-vi[1]*vsum[0]+vi[0]*vsum[1];

			// derivative of center of mass contributions
			// Vector 0 of this face
			Eigen::Vector3d dt = {
				normal[0] * (vsum[0] + vi[0]),
				g[1] * e3[2],
				-g[2] * e3[1]
			};
			dc.col(3 * vinds[0]) += dt;
			dt = {
				-g[0] * e3[2],
				normal[1] * (vsum[1] + vi[1]),
				g[2] * e3[0]
			};
			dc.col(3 * vinds[0] + 1) += dt;
			dt = {
				g[0] * e3[1],
				-g[1] * e3[0],
				normal[2] * (vsum[2] + vi[2])
			};
			dc.col(3 * vinds[0] + 2) += dt;

			// Vector 1 of this face
			dt = {
				normal[0] * (vsum[0] + vj[0]),
				g[1] * e2[2],
				-g[2] * e2[1]
			};
			dc.col(3 * vinds[1]) += dt;
			dt = {
				-g[0] * e2[2],
				normal[1] * (vsum[1] + vj[1]),
				g[2] * e2[0]
			};
			dc.col(3 * vinds[1] + 1) += dt;
			dt = {
				g[0] * e2[1],
				-g[1] * e2[0],
				normal[2] * (vsum[2] + vj[2])
			};
			dc.col(3 * vinds[1] + 2) += dt;

			// Vector 2 of this face
			dt = {
				normal[0] * (vsum[0] + vk[0]),
				g[1] * e1[2],
				-g[2] * e1[1]
			};
			dc.col(3 * vinds[2]) += dt;
			dt = {
				-g[0] * e1[2],
				normal[1] * (vsum[1] + vk[1]),
				g[2] * e1[0]
			};
			dc.col(3 * vinds[2] + 1) += dt;
			dt = {
				g[0] * e1[1],
				-g[1] * e1[0],
				normal[2] * (vsum[2] + vk[2])
			};
			dc.col(3 * vinds[2] + 2) += dt;
		}
		cog_computed = true;
		for (int l = 0;l<dc.cols();l++){
			dcf.col(l) += 3*(dc.col(l)-c*dm[l]/m)/(4*m);
		}
		m /= 18.0; // Authors divide by 6 which is wrong
		c /= (24.0 * m);
		std::cout << "mass: " << m << std::endl;
		std::cout << "COG: " << c << std::endl;
		//std::cout << "Derivative of com sample: " << dcf << std::endl;
		draw_handles();
	}

	bool mouse_down(int button, int modifier) {
		if (button == 0 && create_handle_mode != 0) {
			int v = get_closest_mesh_vertex();
			if (v >= 0) {
				Eigen::RowVector3d pos = viewer->data().V.row(v);
				Eigen::MatrixXd H;
				handles.transfored_handles(H);
				int best = 0;
				for (int i = 0; i < H.rows(); ++i) {
					Eigen::RowVector3d h = H.row(i);
					if ((h - pos).norm() < (Eigen::RowVector3d(H.row(best)) - pos).norm()) {
						best = i;
					}
				}
				moving_handle_id = best;
				sel_pos = H.row(best);
				draw_handles();
				draw_support_polygon();
				return true;
			}
		}
		return false;
	}

	bool mouse_up(int button, int modifier)
	{
		moving_handle_id = -1;
		draw_handles();
		if (support_convex_hull.size() > 2) draw_support_polygon();
		// Check to see if click or rotation
		double dx = viewer->current_mouse_x - viewer->down_mouse_x;
		double dy = viewer->current_mouse_y - viewer->down_mouse_y;
		double dist_sqr = dx * dx + dy * dy;
		if (dist_sqr > 2.0) {
			return false;
		}

		// Left-Click to Add a Handle
		if (button == 0)
		{
			if (viewer->data().V.rows() > 0 && viewer->data().F.rows() > 0 && create_handle_mode == 0) {
				int v = get_closest_mesh_vertex();
				if (v >= 0) {
					Eigen::Vector3d pos = viewer->data().V.row(v);
					handles.add_point_handle(pos);
					draw_handles();
					if (support_convex_hull.size() > 2) draw_support_polygon();
					return true;
				}
			}
		}
		return false;
	}

	bool mouse_move(int mouse_x, int mouse_y) {
		if (create_handle_mode == 1 && moving_handle_id >= 0) {
			// CoG no longer valid so clear it out
			cog_computed = false;

			float x = viewer->current_mouse_x;
			float y = viewer->core().viewport(3) - viewer->current_mouse_y;

			Eigen::RowVector3f orig_pos = sel_pos.cast<float>();

			Eigen::RowVector3f orig_screen_pos;

			igl::project(
				orig_pos,
				viewer->core().view,
				viewer->core().proj,
				viewer->core().viewport,
				orig_screen_pos
			);

			Eigen::RowVector3f new_screen_pos((float)x, (float)y, orig_screen_pos(2));
			Eigen::RowVector3f new_pos;

			igl::unproject(
				new_screen_pos,
				viewer->core().view,
				viewer->core().proj,
				viewer->core().viewport,
				new_pos
			);

			Eigen::RowVector3d pos = new_pos.cast<double>();

			handles.move_handle(
				moving_handle_id,
				pos);
			
			
			viewer->data().set_vertices(lbs_mat * handles.transform());
			draw_handles();
			return true;
		}
		return false;
	}

	int get_closest_mesh_vertex() {
		int index = -1;
		int fid;
		Eigen::Vector3f bc;
		double x = viewer->current_mouse_x;
		double y = viewer->core().viewport(3) - viewer->current_mouse_y;
		Eigen::RowVector3f last_mouse(x, y, 0);
		if (igl::unproject_onto_mesh(
			last_mouse.head(2),
			viewer->core().view,
			viewer->core().proj,
			viewer->core().viewport,
			viewer->data().V,
			viewer->data().F,
			fid,
			bc))
		{
			const Eigen::MatrixXi& F = viewer->data().F;
			int coord;
			bc.maxCoeff(&coord);
			index = F(fid, coord);
		}
		return index;
	}
	
	void set_support_vertices() {
		double minheight = std::numeric_limits<double>::max();
		double maxheight = std::numeric_limits<double>::min();
		for (int i = 0; i < V.rows(); i++) {
			if (V.row(i)[1] < minheight) minheight = V.row(i)[1];
			if (V.row(i)[1] > maxheight) maxheight = V.row(i)[1];
		}
		min_vertex_height = minheight;
		max_vertex_height = maxheight;
		std::cout << "MIN HEIGHT: " << minheight << std::endl;
		
		// Include thrheshold just in case we have other vertices super close to base plane
		double threshold = support_threshold;
		support_vertices.clear();
		for (int i = 0; i < V.rows(); i++) {
			if (V.row(i)[1] < minheight + threshold) support_vertices.push_back(i);
		}

		// Sort this list first by x then y for convex hull calculation
		auto sort_key = [&](int a, int b) {
			if (V.row(a)[0] < V.row(b)[0]) return true;
			if (V.row(a)[0] == V.row(b)[0]) return V.row(a)[2] < V.row(b)[2];
			return false;
		};
		std::sort(support_vertices.begin(), support_vertices.end(), sort_key);

		std::cout << "Support Vertices: \n";
		for (auto i : support_vertices) {
			std::cout << V.row(i) << std::endl;
		}

		compute_support_convex_hull();
	}
	
	// Determine if points a, b, and c make a counter clockwise turn
	bool ccw(int a, int b, int c) {
		auto pa = V.row(a);
		auto pb = V.row(b);
		auto pc = V.row(c);
		
		auto temp = (pb[0] - pa[0]) * (pc[2] - pa[2]) - (pc[0] - pa[0]) * (pb[2] - pa[2]);
		return temp > 0.;
	}

	void draw_handles() {
		Eigen::MatrixXd points, point_colors, line_colors;
		Eigen::MatrixXi lines;
		handles.visualize_handles(points, point_colors, lines, line_colors, create_handle_mode == 1);
		if (moving_handle_id >= 0) {
			point_colors.row(moving_handle_id) = Eigen::Vector3d(1.0, 0.7, 0.2);
		}
		// Draw center of gravity
		if (cog_computed) {
			int oldrows = points.rows();
			int oldcols = points.cols();
			points.conservativeResize(oldrows+1, oldcols);
			points.row(oldrows) = c;
			point_colors.conservativeResize(oldrows + 1, oldcols);
			point_colors.row(oldrows) = Eigen::Vector3d({ 0.0, 0.5, 0.5 });
		}

		viewer->data().clear_labels();
		viewer->data().set_points(points, point_colors);
		viewer->data().set_edges(points, lines, line_colors);
}

	void draw_support_polygon() {
		Eigen::MatrixXd P, P1, P2;
		P1.resize(support_convex_hull.size(), 3);
		P2.resize(support_convex_hull.size(), 3);

		P1.row(0) = V.row(support_convex_hull[0]);

		for (int i = 1; i < support_convex_hull.size(); i++) {
			P1.row(i) = V.row(support_convex_hull[i]);
			P2.row(i - 1) = V.row(support_convex_hull[i]);
		}

		P2.row(support_convex_hull.size() - 1) = P1.row(0);

		for (int i = 0; i < support_convex_hull.size(); i++) {
			P1.row(i)[1] = min_vertex_height;
			P2.row(i)[1] = min_vertex_height;
		}

		viewer->data().add_edges(P1, P2, Eigen::RowVector3d(0., 1., 0.));
	}

	bool will_it_stand(Eigen::Vector3d center_of_gravity) {
		for (int i = 0; i < support_convex_hull.size(); i++) {
			Eigen::Vector3d point = V.row(support_convex_hull[i]);
			Eigen::Vector3d dir_to_support_center = support_center - point;
			Eigen::Vector3d dir_to_center_of_gravity = center_of_gravity - point;
			point[1] = 0.; dir_to_support_center[1] = 0.; dir_to_center_of_gravity[1] = 0.;
			if (dir_to_support_center.dot(dir_to_center_of_gravity) < 0) return false;
		}
		return true;
	}

	void compute_support_convex_hull() {
		// Andrew's monotone chain convex hull algorithm
		if (support_vertices.size() <= 2) return; // Do not recompute unless new support is valid

		std::vector<int> U, L; // indices of vertices in upper and lower hull

		// Lower hull first
		for (int i = 0; i < support_vertices.size(); i++) {
			while (L.size() >= 2 &&
				!ccw(*(L.end() - 2), L.back(), support_vertices[i])) {
				L.pop_back();
			}
			L.push_back(support_vertices[i]);
		}

		// Upper hull
		for (int i = support_vertices.size() - 1; i >= 0; i--) {
			while (U.size() >= 2 &&
				!ccw(*(U.end() - 2), U.back(), support_vertices[i])) {
				U.pop_back();
			}
			U.push_back(support_vertices[i]);
		}

		// Delete repeats
		U.pop_back(); L.pop_back();

		support_convex_hull.clear();
		support_convex_hull = std::vector<int>(L.begin(), L.end());
		support_convex_hull.insert(support_convex_hull.end(), U.begin(), U.end());

		support_center.setZero();
		for (int i = 0; i < support_convex_hull.size(); i++) {
			support_center += V.row(support_convex_hull[i]);
		}
		support_center /= support_convex_hull.size();
	}

	static const int MAX_FACES = 10000;

private:

	igl::opengl::glfw::imgui::ImGuiMenu menu;
    Eigen::Vector3d sc = Eigen::Vector3d({ 1.0, 1.0, 0.0 });//target center of support, for now using 0
	Eigen::MatrixXd V; // Vertex Positions
	Eigen::MatrixXd deformed_V; // Theoretically this is the most up to date version of the vertices after a call to compute_cog
	Eigen::MatrixXi T; // Tetrahedral Elements
	Eigen::MatrixXi F; // Triangular Faces of exterior surface
	float support_threshold = 0.5; // Determine what amount of bottom points are considered support points
	double min_vertex_height = 0., max_vertex_height = 1.;
	std::vector<int> support_vertices; // Indices into V for vertices that will define a support plane
	std::vector<int> support_convex_hull; // Convex hull that represents support polygon
	Eigen::Vector3d support_center; // Center of support polygon
	double m;
	Eigen::Vector3d c;
	Eigen::Matrix<double, 1, Eigen::Dynamic> dm; // (1, 3 * n_vertices) vertex v dm[3*v], dm[3*v + 1], dm[3*v + 2]
	Eigen::Matrix<double, 3, Eigen::Dynamic> dc; // (3, 3 * n_vertices)	
	Eigen::Matrix<double, 3, Eigen::Dynamic> dcf; // (3, 3 * n_vertices)	
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dch; // (3,3*n_handles)
	bool cog_computed = false;
	bool started_mesh_manipulation = false;

	int weight_type = 0; // Type of weights to compute
	                     // 0 = Nearest Neighbor
	                     // 1 = Linear
	                     // 2 = Bounded Biharmonic

	Eigen::MatrixXd lbs_mat; // Linear Blend Skinning Matrix


	// Handle Creation State
	int create_handle_mode = 0; // 0 = add_handle, 1 = manipulate handles
	WeightHandles handles;


	// Handle Manipulation State
	int moving_handle_id = -1;
	Eigen::RowVector3d sel_pos; // Position of handle when it was selected
	//bool fixed = false;
	int n_handles = 0;
};
