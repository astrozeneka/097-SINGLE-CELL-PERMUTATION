import Link from "next/link";

export default function Home() {
  return (
    <div className="min-h-screen bg-slate-950 text-slate-100 p-6">
      <div className="max-w-5xl mx-auto">
        <header className="mb-8">
          <h1 className="text-2xl font-light text-slate-100 tracking-wide">
            Single Cell Analysis Tools
          </h1>
          <div className="h-px bg-gradient-to-r from-slate-700 via-slate-600 to-transparent mt-2"></div>
        </header>

        <div className="space-y-8">
          {/* Neighborhood Analysis Pipeline */}
          <section>
            <h2 className="text-lg font-light text-slate-300 mb-3 tracking-wide">
              Neighborhood Analysis Pipeline
            </h2>
            <div className="grid gap-2">
              <Link
                href="/nhood/001_nhood_run_scimap"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 1: Compute neighborhood
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Define neighborhood using KNN or radius method
                </p>
              </Link>

              <Link
                href="/nhood/002_nhood_cluster"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 2: Cluster motifs
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Identify and cluster neighborhood motifs
                </p>
              </Link>

              <Link
                href="/nhood/003_plot_side_by_side"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 3: Visualize
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Generate side-by-side visualization plots
                </p>
              </Link>
            </div>
          </section>

          {/* Cell-to-Cell Interaction Pipeline */}
          <section>
            <h2 className="text-lg font-light text-slate-300 mb-3 tracking-wide">
              Cell-to-Cell Interaction Pipeline
            </h2>
            <div className="grid gap-2">
              <Link
                href="/permutation-test"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 1: Spatial Distance Permutation Test
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Run permutation tests on spatial distance data
                </p>
              </Link>

              <Link
                href="/spatial-distance"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 2: Spatial Distance Analysis
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Calculate spatial distances between cell phenotypes
                </p>
              </Link>

              <Link
                href="/cell-cell-interaction-analysis"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 3: Cell-Cell Interaction Analysis
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Analyze attraction and avoidance between cell types
                </p>
              </Link>

              <Link
                href="/attraction-avoidance-heatmap"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 4: Attraction-Avoidance Heatmap
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Generate hierarchically clustered heatmaps from permutation test results
                </p>
              </Link>
            </div>
          </section>

          {/* Proximity Analysis Pipeline */}
          <section>
            <h2 className="text-lg font-light text-slate-300 mb-3 tracking-wide">
              Proximity Analysis Pipeline
            </h2>
            <div className="grid gap-2">
              <Link
                href="/proximity-analysis/001_proximity_compute"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 1: Compute
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Compute pairwise proximity scores from spatial coordinates
                </p>
              </Link>

              <Link
                href="/proximity-analysis/002_proximity_analyze"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 2: Analyze
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Compare proximity scores between two groups with statistical testing
                </p>
              </Link>
            </div>
          </section>

          {/* Spatial Score Pipeline */}
          <section>
            <h2 className="text-lg font-light text-slate-300 mb-3 tracking-wide">
              Spatial Score Pipeline
            </h2>
            <div className="grid gap-2">
              <Link
                href="/spatial-score/nearest-cell-distance"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 1: Nearest Cell Distance
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Compute per-cell distances to each phenotype using spatial coordinates
                </p>
              </Link>

              <Link
                href="/spatial-score/compute-ratios"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 2: Compute Ratios
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Compute reference/effector distance ratios per target phenotype and image
                </p>
              </Link>

              <Link
                href="/spatial-score/003_stats"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 3: Stats
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Pairwise Mann-Whitney U test across groups
                </p>
              </Link>

              <Link
                href="/spatial-score/004_boxplot"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Step 4: Boxplot
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Visualize selected target/effector pairs as grouped boxplots
                </p>
              </Link>
            </div>
          </section>

          {/* Other Tools */}
          <section>
            <h2 className="text-lg font-light text-slate-300 mb-3 tracking-wide">
              Other Tools
            </h2>
            <div className="grid gap-2">
              <Link
                href="/attribution-matrix"
                className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
              >
                <h3 className="text-sm font-medium text-slate-100">
                  Attribution Matrix
                </h3>
                <p className="text-xs text-slate-400 mt-0.5">
                  Analyze attribution matrix data
                </p>
              </Link>
            </div>
          </section>
        </div>
      </div>
    </div>
  );
}
