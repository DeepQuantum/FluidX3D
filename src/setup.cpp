#include "setup.hpp"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <vector>

struct Shot {
	string name;
	float3 pos0, pos1, target0, target1;
	float fov;
};

static const vector<Shot> shots = {
	{ "A", float3( 60.0f, -54.0f, 31.0f), float3( 56.0f, -50.0f, 29.0f), float3(0.0f, 11.0f, -1.0f), float3(0.0f, 13.0f, -1.0f), 42.0f },
	{ "B", float3(-34.0f,  78.0f, 14.0f), float3(-34.0f,  78.0f, 14.0f), float3(0.0f,  0.0f,  0.0f), float3(0.0f,  0.0f,  0.0f), 38.0f },
	{ "C", float3( 80.0f,  14.0f,  3.0f), float3( 80.0f,  14.0f,  3.0f), float3(0.0f, 14.0f,  0.0f), float3(0.0f, 14.0f,  0.0f), 52.0f },
	{ "D", float3( 20.0f,  10.0f, 74.0f), float3( 20.0f,  22.0f, 74.0f), float3(0.0f, 14.0f,  0.0f), float3(0.0f, 18.0f,  0.0f), 54.0f },
	{ "E", float3( 44.0f, -20.0f,-30.0f), float3( 44.0f, -20.0f,-30.0f), float3(0.0f, 14.0f,  0.0f), float3(0.0f, 14.0f,  0.0f), 44.0f },
	{ "F", float3( 20.0f, -60.0f, 14.0f), float3( 20.0f, -60.0f, 14.0f), float3(0.0f, 10.0f,  0.0f), float3(0.0f, 10.0f,  0.0f), 34.0f },
	{ "G", float3(  6.0f,  88.0f, 30.0f), float3(  6.0f,  88.0f, 30.0f), float3(0.0f,  0.0f,  0.0f), float3(0.0f,  0.0f,  0.0f), 38.0f },
	{ "H", float3(-44.0f,  26.0f, 62.0f), float3(-44.0f,  26.0f, 62.0f), float3(0.0f, 16.0f,  0.0f), float3(0.0f, 16.0f,  0.0f), 46.0f },
	{ "I", float3( 34.0f,  -6.0f, 16.0f), float3( 34.0f,  -6.0f, 16.0f), float3(11.0f, 26.0f,  0.0f), float3(11.0f, 26.0f,  0.0f), 44.0f },
	{ "J", float3( 20.0f,  92.0f,  7.0f), float3( 20.0f,  92.0f,  7.0f), float3(0.0f,  10.0f,  0.0f), float3(0.0f, 10.0f,  0.0f), 26.0f },
	{ "K", float3( 20.0f,  16.0f, 74.0f), float3( 20.0f,  16.0f, 74.0f), float3(0.0f, 16.0f,  0.0f), float3(0.0f, 16.0f,  0.0f), 54.0f },
	{ "L", float3( 20.0f, -20.0f, 13.0f), float3( 20.0f, -20.0f, 13.0f), float3(7.0f,  2.0f,  0.0f), float3(7.0f,  2.0f,  0.0f), 46.0f },
	{ "M", float3( 26.0f,  38.0f, 14.0f), float3( 26.0f,  38.0f, 14.0f), float3(4.0f, 18.0f,  1.0f), float3(4.0f, 18.0f,  1.0f), 46.0f },
	{ "N", float3( 18.0f, -14.0f,-14.0f), float3( 18.0f, -14.0f,-14.0f), float3(6.0f,  2.0f, -1.0f), float3(6.0f,  2.0f, -1.0f), 46.0f },
	{ "P", float3( 30.0f, -43.0f, 18.0f), float3( 30.0f, -43.0f, 18.0f), float3(1.0f,  1.0f, -1.5f), float3(1.0f,  1.0f, -1.5f), 46.0f },
	{ "Q", float3( 33.3f, -47.8f,-21.2f), float3( 33.3f, -47.8f,-21.2f), float3(1.0f,  1.0f, -1.5f), float3(1.0f,  1.0f, -1.5f), 46.0f }
};

static string env_s(const char* key, const string& fallback) {
	const char* v = std::getenv(key);
	return v ? string(v) : fallback;
}
static float env_f(const char* key, const float fallback) {
	const char* v = std::getenv(key);
	return v ? (float)std::atof(v) : fallback;
}
static ulong env_u(const char* key, const ulong fallback) {
	const char* v = std::getenv(key);
	return v ? (ulong)std::strtoull(v, nullptr, 10) : fallback;
}

static void voxelize_parity(LBM& lbm, const Mesh* mesh, const uchar flag, ulong& solid_cells, ulong& odd_columns, ulong& planform_cells, ulong& frontal_cells) {
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz();
	vector<vector<float>> hits((size_t)Nx*(size_t)Nz);
	const double ox=0.000137, oz=0.000271;
	for(uint i=0u; i<mesh->triangle_number; i++) {
		const double ax=mesh->p0[i].x, ay=mesh->p0[i].y, az=mesh->p0[i].z;
		const double bx=mesh->p1[i].x, by=mesh->p1[i].y, bz=mesh->p1[i].z;
		const double cx=mesh->p2[i].x, cy=mesh->p2[i].y, cz=mesh->p2[i].z;
		const int x0=max((int)std::ceil(std::min({ax, bx, cx})-ox), 0), x1=min((int)std::floor(std::max({ax, bx, cx})-ox), (int)Nx-1);
		const int z0=max((int)std::ceil(std::min({az, bz, cz})-oz), 0), z1=min((int)std::floor(std::max({az, bz, cz})-oz), (int)Nz-1);
		for(int z=z0; z<=z1; z++) {
			for(int x=x0; x<=x1; x++) {
				const double px=(double)x+ox, pz=(double)z+oz;
				const double e0=(bx-ax)*(pz-az)-(bz-az)*(px-ax);
				const double e1=(cx-bx)*(pz-bz)-(cz-bz)*(px-bx);
				const double e2=(ax-cx)*(pz-cz)-(az-cz)*(px-cx);
				if((e0>0.0&&e1>0.0&&e2>0.0)||(e0<0.0&&e1<0.0&&e2<0.0)) {
					hits[(size_t)x+(size_t)z*(size_t)Nx].push_back((float)((e1*ay+e2*by+e0*cy)/(e0+e1+e2)));
				}
			}
		}
	}
	std::atomic<ulong> solid(0ull), odd(0ull), frontal(0ull);
	parallel_for((ulong)Nx*(ulong)Nz, [&](ulong c) {
		vector<float>& h = hits[c];
		if(h.empty()) return;
		frontal++;
		std::sort(h.begin(), h.end());
		if(h.size()%2u) odd++;
		const uint x=(uint)(c%(ulong)Nx), z=(uint)(c/(ulong)Nx);
		ulong count = 0ull;
		for(size_t k=0u; k+1u<h.size(); k+=2u) {
			const int y0=max((int)std::ceil(h[k]), 0), y1=min((int)std::floor(h[k+1u]), (int)Ny-1);
			for(int y=y0; y<=y1; y++) {
				lbm.flags[lbm.index(x, (uint)y, z)] = flag;
				count++;
			}
		}
		solid += count;
	});
	vector<uchar> plan((size_t)Nx*(size_t)Ny, (uchar)0u);
	for(ulong n=0ull; n<lbm.get_N(); n++) {
		if(lbm.flags[n]&TYPE_S) {
			uint x=0u, y=0u, z=0u;
			lbm.coordinates(n, x, y, z);
			plan[(size_t)x+(size_t)y*(size_t)Nx] = (uchar)1u;
		}
	}
	ulong planform = 0ull;
	for(const uchar p : plan) planform += p;
	solid_cells=solid.load(); odd_columns=odd.load(); planform_cells=planform; frontal_cells=frontal.load();
}

void main_setup() {
	const string stl = env_s("VS_STL", "");
	const string out = env_s("VS_OUT", "");
	if(stl.empty()||out.empty()) print_error("VS_STL and VS_OUT must be set.");
	const float aoa = env_f("VS_AOA", 5.0f);
	const uint memory = (uint)env_u("VS_MB", 400ull);
	const float lbm_u = env_f("VS_U", 0.075f);
	const float re = env_f("VS_RE", 1.0E7f);
	const float si_u = env_f("VS_SI_U", 77.0f);
	const float si_rho = 1.225f;
	const float si_length = env_f("VS_SI_LENGTH", 44.8f);
	const float box_x=env_f("VS_BOX_X", 48.0f), box_y=env_f("VS_BOX_Y", 120.0f), box_z=env_f("VS_BOX_Z", 32.0f);
	const float nose_gap = env_f("VS_NOSE_GAP", 14.0f);
	const float z_frac = env_f("VS_Z_FRAC", 0.42f);
	const ulong lbm_T = env_u("VS_T", 2000ull);
	const ulong frame_every = max(env_u("VS_FRAME_EVERY", 500ull), 1ull);
	const ulong first_frame = env_u("VS_FIRST_FRAME", 0ull);
	const ulong force_every = max(env_u("VS_FORCE_EVERY", 200ull), frame_every);
	const ulong check_every = env_u("VS_CHECK_EVERY", 2000ull);
	const ulong gpu_check_every = env_u("VS_GPU_CHECK_EVERY", 0ull);
	if(env_u("VS_DROP_HOST_FIELDS", 0ull)&&check_every) print_error("Host-field release requires VS_CHECK_EVERY=0.");
	const string shot_names = env_s("VS_SHOTS", "A");
	const string bare_names = env_s("VS_BARE", "");
	const string sweep_names = env_s("VS_SWEEP_SHOTS", "");
	const string q_sweep = env_s("VS_Q_SWEEP", "");
	const ulong max_frames = env_u("VS_MAX_FRAMES", 6000ull);
	const int film_modes = (int)env_u("VS_VIS", (ulong)(VIS_FLAG_SURFACE|VIS_Q_CRITERION));

	const uint domains_y = (uint)env_u("VS_DY", 1ull);
	const float target_dx = env_f("VS_DX", 0.0f);
	uint3 lbm_N = target_dx>0.0f ? uint3(to_uint(box_x/target_dx), to_uint(box_y/target_dx), to_uint(box_z/target_dx)) : resolution(float3(box_x, box_y, box_z), memory);
	if(target_dx>0.0f) lbm_N.y += (domains_y-lbm_N.y%domains_y)%domains_y;
	else lbm_N.y -= lbm_N.y%domains_y;
	const float cells_per_m = (float)lbm_N.x/box_x;
	const float lbm_length = si_length*cells_per_m;
	const float lbm_nu = lbm_u*lbm_length/re;
	units.set_m_kg_s(lbm_length, lbm_u, 1.0f, si_length, si_u, si_rho);
	const float si_dt = si_length/lbm_length*lbm_u/si_u;
	print_info("grid "+to_string(lbm_N.x)+"x"+to_string(lbm_N.y)+"x"+to_string(lbm_N.z)+" = "+to_string(to_uint((float)lbm_N.x*(float)lbm_N.y*(float)lbm_N.z*1E-6f))+" Mcells, dx = "+to_string(1.0f/cells_per_m, 4u)+" m, hull length "+to_string(lbm_length, 1u)+" cells, dt = "+to_string(si_dt*1E6f, 2u)+" us");
	print_info("lbm_u = "+to_string(lbm_u, 4u)+", lbm_nu = "+to_string(lbm_nu*1E6f, 4u)+"e-6, Re_L nominal = "+to_string(re, 0u)+", real air at this speed would be "+to_string(si_u*si_length/1.48E-5f, 0u));
	const float rake_x = env_f("VS_RAKE_X", 0.0f);
	if(rake_x>0.0f) {
		const float zc = (z_frac-0.5f)*(float)lbm_N.z;
		const string rx=std::to_string(rake_x*cells_per_m), r0=std::to_string(zc+env_f("VS_RAKE_Z0", -7.5f)*cells_per_m), r1=std::to_string(zc+env_f("VS_RAKE_Z1", 5.5f)*cells_per_m);
		setenv("VS_RAKE_XC", rx.c_str(), 1);
		setenv("VS_RAKE_Z0C", r0.c_str(), 1);
		setenv("VS_RAKE_Z1C", r1.c_str(), 1);
		setenv("VS_SL_DT0", "1.0", 1);
	}
	const bool render_only = env_u("VS_RENDER_ONLY", 0ull)!=0ull;
	LBM lbm(render_only ? uint3(8u, 8u, 8u) : lbm_N, 1u, render_only ? 1u : domains_y, 1u, lbm_nu);
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz();

	const bool tilt_flow = env_u("VS_TILT", 1ull)!=0ull;
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(tilt_flow ? 0.0f : -aoa))*float3x3(float3(0, 0, 1), radians(180.0f));
	const float wind_y = tilt_flow ? cos(radians(aoa)) : 1.0f;
	const float wind_z = tilt_flow ? sin(radians(aoa)) : 0.0f;
	Mesh* mesh = read_stl(stl, cells_per_m, rotation, float3(0.0f));
	mesh->find_bounds();
	const float3 shift = float3(lbm.center().x-0.5f*(mesh->pmin.x+mesh->pmax.x), nose_gap*cells_per_m-mesh->pmin.y, z_frac*(float)Nz-0.5f*(mesh->pmin.z+mesh->pmax.z));
	mesh->translate(shift);
	mesh->find_bounds();
	const float3 ship_center = 0.5f*(mesh->pmin+mesh->pmax);
	print_info("mesh bbox cells min "+to_string(mesh->pmin.x, 1u)+" "+to_string(mesh->pmin.y, 1u)+" "+to_string(mesh->pmin.z, 1u)+" max "+to_string(mesh->pmax.x, 1u)+" "+to_string(mesh->pmax.y, 1u)+" "+to_string(mesh->pmax.z, 1u));
	const float3 box_center = lbm.center();
	if(render_only) {
		lbm.run(0u, 0u);
		lbm.graphics.set_mesh(mesh);
		lbm.graphics.visualization_modes = VIS_MESH;
		for(const Shot& s : shots) {
			if(shot_names.find(s.name)==string::npos) continue;
			const float3 p = ship_center-box_center+cells_per_m*s.pos0;
			const float3 c = ship_center-box_center+cells_per_m*s.target0;
			const float3 d = p-c;
			lbm.graphics.set_camera_free(p, degrees(atan2(d.y, d.x)), degrees(asin(d.z/length(d))), s.fov);
			lbm.graphics.write_frame(out+"bare/"+s.name+"/");
		}
		return;
	}

	Clock vox_clock;
	vox_clock.start();
	ulong solid_cells=0ull, odd_columns=0ull, planform_cells=0ull, frontal_cells=0ull;
	voxelize_parity(lbm, mesh, (uchar)(TYPE_S|TYPE_X), solid_cells, odd_columns, planform_cells, frontal_cells);
	const float si_planform = (float)planform_cells/sq(cells_per_m);
	const float si_frontal = (float)frontal_cells/sq(cells_per_m);
	print_info("voxelized in "+to_string(vox_clock.stop(), 2u)+" s: solid "+to_string(solid_cells)+", odd-parity columns "+to_string(odd_columns)+", frontal "+to_string(si_frontal, 1u)+" m2, planform "+to_string(si_planform, 1u)+" m2, solid volume "+to_string((float)solid_cells/(cells_per_m*cells_per_m*cells_per_m), 1u)+" m3");

	parallel_for(lbm.get_N(), [&](ulong n) {
		uint x=0u, y=0u, z=0u;
		lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]!=(TYPE_S|TYPE_X)) {
			lbm.u.y[n] = lbm_u*wind_y;
			lbm.u.z[n] = lbm_u*wind_z;
		}
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E;
	});
	vector<uchar> host_flags;
	if(check_every) {
		host_flags.resize((size_t)lbm.get_N());
		for(ulong n=0ull; n<lbm.get_N(); n++) host_flags[(size_t)n] = lbm.flags[n];
	}

	lbm.graphics.visualization_modes = film_modes;
	lbm.graphics.field_mode = (int)env_u("VS_FIELD_MODE", 0ull);
	lbm.graphics.slice_mode = (int)env_u("VS_SLICE_MODE", 0ull);
	lbm.graphics.slice_x = (int)((float)Nx*env_f("VS_SLICE_X", 0.5f));
	lbm.graphics.slice_y = (int)((float)Ny*env_f("VS_SLICE_Y", 0.75f));
	lbm.graphics.slice_z = (int)((float)Nz*env_f("VS_SLICE_Z", 0.5f));

	vector<const Shot*> active;
	for(const Shot& s : shots) if(shot_names.find(s.name)!=string::npos) active.push_back(&s);
	const ulong planned = (lbm_T-min(first_frame, lbm_T))/frame_every+1ull;
	print_info("frames planned per shot "+to_string(planned)+", shots "+to_string((uint)active.size()));
	if(planned*(ulong)active.size()>max_frames) print_error("Frame plan exceeds VS_MAX_FRAMES.");

	const auto aim = [&](const Shot& s, const float eased) {
		const float3 p_m = s.pos0+eased*(s.pos1-s.pos0);
		const float3 t_m = s.target0+eased*(s.target1-s.target0);
		const float3 p = ship_center-box_center+cells_per_m*p_m;
		const float3 c = ship_center-box_center+cells_per_m*t_m;
		const float3 d = p-c;
		lbm.graphics.set_camera_free(p, degrees(atan2(d.y, d.x)), degrees(asin(d.z/length(d))), s.fov);
	};
	const auto stills = [&](const int modes, const string& names, const string& folder, const float eased) {
		lbm.graphics.visualization_modes = modes;
		for(const Shot& s : shots) {
			if(names.find(s.name)==string::npos) continue;
			aim(s, eased);
			lbm.graphics.write_frame(out+folder+"/"+s.name+"/");
		}
		lbm.graphics.visualization_modes = film_modes;
	};

	create_folder(out);
#ifdef FORCE_FIELD
	const string log_path = out+"forces.csv";
	FILE* force_log = fopen(log_path.c_str(), "w");
	if(force_log==nullptr) print_error("Cannot open force log.");
	fprintf(force_log, "t,si_time_s,Fx_N,Fy_drag_N,Fz_lift_N,CD_planform,CL_planform,Mroll_Nm,Mpitch_Nm,Myaw_Nm\n");
	const float3 mrp_mesh = cells_per_m*(rotation*float3(env_f("VS_MRP_X", 0.0f), env_f("VS_MRP_Y", 0.0f), env_f("VS_MRP_Z", 0.0f)))+shift;
	const float3 mrp_box = mrp_mesh+float3(0.5f-0.5f*(float)Nx, 0.5f-0.5f*(float)Ny, 0.5f-0.5f*(float)Nz);
	if(lbm.get_D()>1u) print_warning("Moments are only valid for a single domain.");
	print_info("planform "+to_string(si_planform, 2u)+" m2, frontal "+to_string(si_frontal, 2u)+" m2, moment reference box cells "+to_string(mrp_box.x, 1u)+" "+to_string(mrp_box.y, 1u)+" "+to_string(mrp_box.z, 1u));
	const float q_ref = 0.5f*si_rho*sq(si_u)*si_planform;
#endif

	lbm.run(0u, lbm_T);
	lbm.graphics.set_mesh(mesh);
	if(!bare_names.empty()) {
		stills(VIS_FLAG_SURFACE, bare_names, "bare_vox", 0.0f);
		stills(VIS_MESH, bare_names, "bare_mesh", 0.0f);
	}
	bool blown = false;
	ulong next_gpu_check = 0ull;
	Clock wall;
	wall.start();
	double rate_wall = 0.0;
	ulong rate_t=lbm.get_t(), rate_up=lbm.out_of_core_uploaded(), rate_down=lbm.out_of_core_downloaded();
	while(!blown) {
		const ulong t = lbm.get_t();
		if(gpu_check_every&&(t>=next_gpu_check||t>=lbm_T)) {
			lbm.update_fields();
			float r_min=2.0f, r_max=0.0f, u_max=0.0f;
			uint bad=0u;
			for(uint d=0u; d<lbm.get_D(); d++) {
				const auto state = lbm.lbm_domain[d]->sampled_state();
				r_min=fmin(r_min, state.rho_min); r_max=fmax(r_max, state.rho_max); u_max=fmax(u_max, state.speed_max); bad+=state.bad;
			}
			fprintf(stdout, "\nsampled check t=%llu rho %.5f..%.5f |u|max %.4f bad_samples %u wall %.1f s\n", (unsigned long long)t, r_min, r_max, u_max, bad, wall.stop());
			fflush(stdout);
			if(bad||u_max>0.5f) blown = true;
			if(lbm.out_of_core()&&t>rate_t) {
				const double now = wall.stop(), seconds = now-rate_wall, steps = (double)(t-rate_t);
				const double up = (double)(lbm.out_of_core_uploaded()-rate_up), down = (double)(lbm.out_of_core_downloaded()-rate_down);
				fprintf(stdout, "out-of-core window t=%llu..%llu: %.3f s/step, %.0f MLUPs, upload %.1f GB/s, download %.1f GB/s (frames and checks included)\n", (unsigned long long)rate_t, (unsigned long long)t, seconds/steps, 1E-6*(double)lbm.get_N()*steps/seconds, 1E-9*up/seconds, 1E-9*down/seconds);
				fflush(stdout);
				rate_wall=now; rate_t=t; rate_up=lbm.out_of_core_uploaded(); rate_down=lbm.out_of_core_downloaded();
			}
			next_gpu_check = t+(t==0ull ? min(gpu_check_every, 5ull*frame_every) : gpu_check_every);
		}
#ifdef FORCE_FIELD
		if(t%force_every==0ull) {
			const float3 f = lbm.object_force((uchar)(TYPE_S|TYPE_X));
			const float3 m = lbm.object_torque(mrp_box, (uchar)(TYPE_S|TYPE_X));
			const float fx=units.si_F(f.x), fy=units.si_F(f.y*wind_y+f.z*wind_z), fz=units.si_F(f.z*wind_y-f.y*wind_z);
			const float mx=-units.si_F(m.x)/cells_per_m, my=-units.si_F(m.y)/cells_per_m, mz=units.si_F(m.z)/cells_per_m;
			fprintf(force_log, "%llu,%.5f,%.1f,%.1f,%.1f,%.4f,%.4f,%.1f,%.1f,%.1f\n", (unsigned long long)t, (double)t*(double)si_dt, fx, fy, fz, fy/q_ref, fz/q_ref, my, mx, mz);
			fflush(force_log);
			if(!(fy==fy)||!(fz==fz)) blown = true;
		}
#endif
		if(t%frame_every==0ull&&t>=first_frame) {
			const float progress = lbm_T>first_frame ? (float)(t-first_frame)/(float)(lbm_T-first_frame) : 0.0f;
			const float eased = progress*progress*(3.0f-2.0f*progress);
			for(const Shot* s : active) {
				aim(*s, eased);
				lbm.graphics.write_frame(out+s->name+"/");
			}
		}
		if(check_every>0ull&&t%check_every==0ull) {
			lbm.update_fields();
			lbm.rho.read_from_device();
			lbm.u.read_from_device();
			float rho_min=2.0f, rho_max=0.0f, u_max=0.0f;
			ulong bad = 0ull;
			for(ulong n=0ull; n<lbm.get_N(); n++) {
				if(host_flags[(size_t)n]&TYPE_S) continue;
				const float r=lbm.rho[n], ux=lbm.u.x[n], uy=lbm.u.y[n], uz=lbm.u.z[n];
				const float u2=ux*ux+uy*uy+uz*uz;
				if(!(r==r)||!(u2==u2)) { bad++; continue; }
				rho_min=fmin(rho_min, r); rho_max=fmax(rho_max, r); u_max=fmax(u_max, u2);
			}
			fprintf(stdout, "\ncheck t=%llu rho %.5f..%.5f |u|max %.4f nan_cells %llu wall %.1f s\n", (unsigned long long)t, rho_min, rho_max, sqrt(u_max), (unsigned long long)bad, wall.stop());
			fflush(stdout);
			if(bad>0ull||sqrt(u_max)>0.5f) blown = true;
		}
		if(t>=lbm_T) break;
		lbm.run(min(frame_every, lbm_T-t), lbm_T);
	}
	if(env_u("VS_HASH_FIELDS", 0ull)&&!blown) lbm.print_field_hashes();
	if(!q_sweep.empty()&&!blown) {
		const float q_keep = vs_q_runtime;
		std::stringstream list(q_sweep);
		string item;
		while(std::getline(list, item, ',')) {
			vs_q_runtime = (float)std::atof(item.c_str());
			stills(VIS_Q_CRITERION|(film_modes&(VIS_FLAG_SURFACE|VIS_MESH)), sweep_names, "q"+item, 1.0f);
		}
		vs_q_runtime = q_keep;
	}
#ifdef FORCE_FIELD
	fclose(force_log);
#endif
	if(blown) println("\nBLOWUP detected, run stopped.");
	println("\nfinished at t="+to_string(lbm.get_t())+" wall "+to_string(wall.stop(), 1u)+" s");
}
