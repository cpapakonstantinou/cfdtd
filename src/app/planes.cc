/**
 * @file app/planes.cc
 * @brief Where the boxes and planes of a run sit in the global lattice, app/planes.h
 *
 * Application policy: the TF/SF box, ART planes and flux boxes from the [boxes] / analysis settings,
 * and the plot planes, some of which follow the antenna shape. The grid module maps the result onto
 * the ranks, grid::set_lattice_data. Moved out of the decomposition on 2026-11-13, same arithmetic.
 */
#include "app/planes.h"
#include "app/run_config.h"
#include "geometry/shapes.h"

fdtd::grid::GlobalPlanes fdtd::app::global_planes(const fdtd::grid::Lattice &lat)
{
	const RunConfig &run_cfg = run_config();
	grid::GlobalPlanes p;
	p.tfsf_nea[0] = lat.cpml[0] + lat.tfsf_inset[0];	// it could have been called TFSFglobneaX
	p.tfsf_nea[1] = lat.cpml[1] + lat.tfsf_inset[1];
	p.tfsf_nea[2] = lat.cpml[2] + lat.tfsf_inset[2];
	p.tfsf_far[0] = lat.cells[0] - (lat.cpml[0] + lat.tfsf_inset[0]);
	p.tfsf_far[1] = lat.cells[1] - (lat.cpml[1] + lat.tfsf_inset[1] + run_cfg.source.tfsf_y_shrink);
	p.tfsf_far[2] = lat.cells[2] - (lat.cpml[2] + lat.tfsf_inset[2]);

	p.art_nea_y = int(lat.cpml[1] + run_cfg.analysis.nr_qsca_y);
	p.art_far_y = int(lat.cells[1] - (lat.cpml[1] + run_cfg.analysis.nr_qtot_y));

	p.qsca_nea[0] = int(lat.cpml[0] + run_cfg.analysis.nr_qsca_x);
	p.qsca_far[0] = int(lat.cells[0] - (lat.cpml[0] + run_cfg.analysis.nr_qsca_x));
	p.qtot_nea[0] = int(lat.cpml[0] + run_cfg.analysis.nr_qtot_x);
	p.qtot_far[0] = int(lat.cells[0] - (lat.cpml[0] + run_cfg.analysis.nr_qtot_x));

	p.qsca_nea[1] = int(lat.cpml[1] + run_cfg.analysis.nr_qsca_y);
	p.qsca_far[1] = int(lat.cells[1] - (lat.cpml[1] + run_cfg.analysis.nr_qsca_y));
	p.qtot_nea[1] = int(lat.cpml[1] + run_cfg.analysis.nr_qtot_y);
	p.qtot_far[1] = int(lat.cells[1] - (lat.cpml[1] + run_cfg.analysis.nr_qtot_y));

	p.qsca_nea[2] = int(lat.cpml[2] + run_cfg.analysis.nr_qsca_z);
	p.qsca_far[2] = int(lat.cells[2] - (lat.cpml[2] + run_cfg.analysis.nr_qsca_z));
	p.qtot_nea[2] = int(lat.cpml[2] + run_cfg.analysis.nr_qtot_z);
	p.qtot_far[2] = int(lat.cells[2] - (lat.cpml[2] + run_cfg.analysis.nr_qtot_z));

	p.plot_x = int(lat.cpml[0] + run_cfg.geometry.center_position_x);
	p.plot_z = int(lat.cpml[2] + run_cfg.geometry.center_position_z);
	p.plot_y2 = int(lat.cells[1] - lat.cpml[1] - lat.tfsf_inset[1] - 5);
	p.plot_y1 = int(lat.cpml[1] + lat.tfsf_inset[1] / 2);

	if (run_cfg.geometry.antenna_type == 2300) // JOH  cylindric dipole -> Z shift by number of cells of length / 2 
	{
		p.plot_z += int(run_cfg.geometry.length / lat.d[2] / 2);
	}

	if (run_cfg.geometry.antenna_type == fdtd::geometry::shape::dipole || run_cfg.geometry.antenna_type == fdtd::geometry::shape::l_antenna || run_cfg.geometry.antenna_type == fdtd::geometry::shape::lateral_dipole || run_cfg.geometry.antenna_type == fdtd::geometry::shape::round_butterfly || run_cfg.geometry.antenna_type == fdtd::geometry::shape::round_dipole || run_cfg.geometry.antenna_type == fdtd::geometry::shape::one_side_round_dipole || run_cfg.geometry.antenna_type == fdtd::geometry::shape::round_tripole || run_cfg.geometry.antenna_type == fdtd::geometry::shape::round_hybrid || run_cfg.geometry.antenna_type == fdtd::geometry::shape::x_shape || run_cfg.geometry.antenna_type == fdtd::geometry::shape::rounded_monopole || run_cfg.geometry.antenna_type == fdtd::geometry::shape::rounded_dipole || run_cfg.geometry.antenna_type == fdtd::geometry::shape::rounded_dipole_offset_planes)
	{
		p.plot_y = int(lat.cpml[1] + run_cfg.geometry.center_position_y) - int(run_cfg.geometry.thickness / lat.d[1] / 2) - int(run_cfg.geometry.antenna_shift / lat.d[1]);
		if (run_cfg.geometry.graphene_optim != 0)
			p.plot_y = int(lat.cpml[1] + run_cfg.geometry.center_position_y) - int(run_cfg.geometry.antenna_shift / lat.d[1]); // + 1;
	}
	else if (run_cfg.geometry.antenna_type == fdtd::geometry::shape::hole)
		p.plot_y = int(lat.cpml[1] + run_cfg.geometry.center_position_y) - int(0.5 * run_cfg.geometry.middle_layer / lat.d[1]);
	else if (run_cfg.geometry.antenna_type == fdtd::geometry::shape::sphere || run_cfg.geometry.antenna_type == fdtd::geometry::shape::bowtie)
		p.plot_y = int(lat.cpml[1] + run_cfg.geometry.center_position_y - run_cfg.geometry.antenna_shift / lat.d[1]);
	else if (run_cfg.geometry.antenna_type == fdtd::geometry::shape::sphere_dimer || run_cfg.geometry.antenna_type == fdtd::geometry::shape::sphere_trimer || run_cfg.geometry.antenna_type == fdtd::geometry::shape::sphere_quadrumer)
		p.plot_y = int(lat.cpml[1] + run_cfg.geometry.center_position_y - run_cfg.geometry.antenna_shift / lat.d[1]);
	else
		p.plot_y = int(lat.cpml[1] + run_cfg.geometry.center_position_y);

	if (run_cfg.geometry.antenna_type == fdtd::geometry::shape::rounded_monopole || run_cfg.geometry.antenna_type == fdtd::geometry::shape::rounded_dipole_offset_planes || (run_cfg.geometry.antenna_type == fdtd::geometry::shape::rounded_triangle && !run_cfg.materials.hydrodynamics))
	{
		p.plot_y1 = int(lat.cpml[1] + run_cfg.geometry.center_position_y - 3);
		p.plot_y2 = int(lat.cpml[1] + run_cfg.geometry.center_position_y - 5);

	}

	if (run_cfg.geometry.antenna_type == fdtd::geometry::shape::rounded_dipole)
	{
		p.plot_y1 = int(lat.cpml[1] + run_cfg.geometry.center_position_y - 3);
		p.plot_y2 = int(lat.cpml[1] + run_cfg.geometry.center_position_y - 5);

	}

	if (run_cfg.source.mode_source)
	{
		p.plot_y2 = int(lat.cpml[1] + lat.tfsf_inset[1]);
	}

	if (run_cfg.geometry.antenna_type == fdtd::geometry::shape::dipole)
	{
		p.plot_y = int(lat.cpml[1] + run_cfg.geometry.center_position_y - 1);
		p.plot_y1 = int(lat.cpml[1] + run_cfg.geometry.center_position_y - 3);
		p.plot_y2 = int(lat.cpml[1] + run_cfg.geometry.center_position_y - 5);
	}

	if (run_cfg.geometry.antenna_type == fdtd::geometry::shape::round_hybrid)
		p.plot_y1 = int(lat.cpml[1] + run_cfg.geometry.center_position_y - int(0.5 * (run_cfg.geometry.middle_layer) / lat.d[1])); //film-coupled hybrid nanoantenna

	p.plot_y_not_on_rank_start = run_cfg.geometry.graphene_optim;
	return p;
}
