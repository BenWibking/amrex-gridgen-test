//==============================================================================
// TwoMomentRad - a radiation transport library for patch-based AMR codes
// Copyright 2024 Benjamin Wibking.
// Released under the MIT license. See LICENSE file included in the GitHub repo.
//==============================================================================
/// \file testDiskGalaxyGrid.cpp
/// \brief Minimal grid-generation test for the DiskGalaxy refinement geometry.
///

#include <cmath>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <string>

#include "AMReX.H"
#include "AMReX_AmrCore.H"
#include "AMReX_GpuDevice.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"
#include "AMReX_TagBox.H"

namespace
{

class BoxNDDynamic
{
      public:
	friend std::ostream &operator<<(std::ostream &os, BoxNDDynamic const &b);
	BoxNDDynamic(amrex::Box const &b, int dim) : m_box(b), m_dim(dim) {}

      private:
	amrex::Box m_box;
	int m_dim;
};

std::ostream &operator<<(std::ostream &os, BoxNDDynamic const &b)
{
	if (b.m_dim == 1) {
		os << "("
		   << "(" << b.m_box.smallEnd(0) << ")" << " "
		   << "(" << b.m_box.bigEnd(0) << ")" << " "
		   << "(" << b.m_box.type(0) << ")"
		   << ")";
	} else if (b.m_dim == 2) {
		os << "("
		   << "(" << b.m_box.smallEnd(0) << "," << b.m_box.smallEnd(1) << ")" << " "
		   << "(" << b.m_box.bigEnd(0) << "," << b.m_box.bigEnd(1) << ")" << " "
		   << "(" << b.m_box.type(0) << "," << b.m_box.type(1) << ")"
		   << ")";
	} else {
		os << b.m_box;
	}
	return os;
}

} // namespace

class DiskGalaxyGrid final : public amrex::AmrCore
{
      public:
	constexpr static amrex::Real parsec = 3.085677587679311e18; // cm

	DiskGalaxyGrid() = default;

	void InitData() { InitFromScratch(0.0); }

	void PrintGrids() const
	{
		std::string gridfile = "DiskGalaxyGrid.grids";
		amrex::ParmParse pp("disk_galaxy_grid");
		pp.query("gridfile", gridfile);

		if (amrex::ParallelDescriptor::IOProcessor()) {
			std::ofstream ofs(gridfile);
			if (!ofs) {
				amrex::Abort("Failed to open DiskGalaxyGrid gridfile output.");
			}

			ofs << " " << std::setw(2) << finestLevel() + 1 << "\n";
			for (int lev = 0; lev <= finestLevel(); ++lev) {
				amrex::Box const prob_domain = Geom(lev).Domain();
				auto const &ba = boxArray(lev);
				ofs << "   " << BoxNDDynamic(prob_domain, AMREX_SPACEDIM) << "  " << ba.size() << "\n";
				for (int ibox = 0; ibox < static_cast<int>(ba.size()); ++ibox) {
					ofs << "      " << BoxNDDynamic(ba[ibox], AMREX_SPACEDIM) << "\n";
				}
			}
		}
	}

      private:
	void MakeNewLevelFromScratch(int lev, amrex::Real time, amrex::BoxArray const &ba, amrex::DistributionMapping const &dm) override
	{
		amrex::ignore_unused(lev, time, ba, dm);
	}

	void MakeNewLevelFromCoarse(int lev, amrex::Real time, amrex::BoxArray const &ba, amrex::DistributionMapping const &dm) override
	{
		amrex::ignore_unused(lev, time, ba, dm);
	}

	void RemakeLevel(int lev, amrex::Real time, amrex::BoxArray const &ba, amrex::DistributionMapping const &dm) override
	{
		amrex::ignore_unused(lev, time, ba, dm);
	}

	void ClearLevel(int lev) override { amrex::ignore_unused(lev); }

	void ErrorEst(int lev, amrex::TagBoxArray &tags, amrex::Real /*time*/, int /*ngrow*/) override
	{
		amrex::ParmParse const pp("agora_galaxy");
		amrex::Real refine_Rmax_kpc = NAN;
		amrex::Real refine_zmax_kpc = NAN;
		pp.query("refine_Rmax_kpc", refine_Rmax_kpc);
		pp.query("refine_zmax_kpc", refine_zmax_kpc);
		amrex::Real const refine_Rmax = refine_Rmax_kpc * (1.0e3 * parsec);
		amrex::Real const refine_zmax = refine_zmax_kpc * (1.0e3 * parsec);

		auto const prob_lo = Geom(lev).ProbLoArray();
		auto const dx = Geom(lev).CellSizeArray();
		auto const tag = tags.arrays();

		amrex::ParallelFor(tags, [=] AMREX_GPU_DEVICE(int bx, int i, int j, int k) noexcept {
			amrex::Real const x0 = prob_lo[0] + (i * dx[0]);
			amrex::Real const y0 = prob_lo[1] + (j * dx[1]);
			amrex::Real const z0 = prob_lo[2] + (k * dx[2]);

			amrex::Real const x1 = prob_lo[0] + ((i + 1) * dx[0]);
			amrex::Real const y1 = prob_lo[1] + ((j + 1) * dx[1]);
			amrex::Real const z1 = prob_lo[2] + ((k + 1) * dx[2]);

			auto tagIfPointInRegion = [=](amrex::Real x, amrex::Real y, amrex::Real z) {
				amrex::Real const R = std::sqrt(x * x + y * y);
				if ((R < refine_Rmax) && (std::abs(z) < refine_zmax)) {
					tag[bx](i, j, k) = amrex::TagBox::SET;
				}
			};

			for (auto const &x : {x0, x1}) {
				for (auto const &y : {y0, y1}) {
					for (auto const &z : {z0, z1}) {
						tagIfPointInRegion(x, y, z);
					}
				}
			}
		});
		amrex::Gpu::streamSynchronize();
	}
};

auto main(int argc, char* argv[]) -> int
{
  amrex::Initialize(argc, argv);
  {
    DiskGalaxyGrid amr;
    amr.InitData();
    amr.PrintGrids();
  }
  amrex::Finalize();
  return 0;
}
