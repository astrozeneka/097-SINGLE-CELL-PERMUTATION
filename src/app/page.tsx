import Link from "next/link";

export default function Home() {
  return (
    <div className="min-h-screen bg-slate-950 text-slate-100 p-6">
      <div className="max-w-5xl mx-auto">
        <header className="mb-6">
          <h1 className="text-2xl font-light text-slate-100 tracking-wide">
            Single Cell Analysis Tools
          </h1>
          <div className="h-px bg-gradient-to-r from-slate-700 via-slate-600 to-transparent mt-2"></div>
        </header>

        <div className="grid gap-2">
          <Link
            href="/permutation-test"
            className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
          >
            <h2 className="text-sm font-medium text-slate-100">
              Spatial Distance Permutation Test
            </h2>
            <p className="text-xs text-slate-400 mt-0.5">
              Run permutation tests on spatial distance data
            </p>
          </Link>

          <Link
            href="/spatial-distance"
            className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
          >
            <h2 className="text-sm font-medium text-slate-100">
              Spatial Distance Analysis
            </h2>
            <p className="text-xs text-slate-400 mt-0.5">
              Calculate spatial distances between cell phenotypes
            </p>
          </Link>

          <Link
            href="/cell-cell-interaction-analysis"
            className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
          >
            <h2 className="text-sm font-medium text-slate-100">
              Cell-Cell Interaction Analysis
            </h2>
            <p className="text-xs text-slate-400 mt-0.5">
              Analyze attraction and avoidance between cell types
            </p>
          </Link>

          <Link
            href="/attribution-matrix"
            className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
          >
            <h2 className="text-sm font-medium text-slate-100">
              Attribution Matrix
            </h2>
            <p className="text-xs text-slate-400 mt-0.5">
              Analyze attribution matrix data
            </p>
          </Link>

          <Link
            href="/attraction-avoidance-heatmap"
            className="block px-4 py-2.5 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
          >
            <h2 className="text-sm font-medium text-slate-100">
              Attraction-Avoidance Heatmap
            </h2>
            <p className="text-xs text-slate-400 mt-0.5">
              Generate hierarchically clustered heatmaps from permutation test results
            </p>
          </Link>
        </div>
      </div>
    </div>
  );
}