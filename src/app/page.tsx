import Link from "next/link";

export default function Home() {
  return (
    <div className="min-h-screen bg-slate-950 text-slate-100 p-6">
      <div className="max-w-5xl mx-auto">
        <header className="mb-8">
          <h1 className="text-3xl font-light text-slate-100 tracking-wide">
            Single Cell Analysis Tools
          </h1>
          <div className="h-px bg-gradient-to-r from-slate-700 via-slate-600 to-transparent mt-2"></div>
        </header>

        <div className="grid gap-4">
          <Link
            href="/permutation-test"
            className="block p-6 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
          >
            <h2 className="text-xl font-medium text-slate-100 mb-2">
              Spatial Distance Permutation Test
            </h2>
            <p className="text-sm text-slate-400">
              Run permutation tests on spatial distance data
            </p>
          </Link>

          <Link
            href="/attribution-matrix"
            className="block p-6 bg-slate-900 hover:bg-slate-800 border border-slate-800 hover:border-slate-700 transition-colors"
          >
            <h2 className="text-xl font-medium text-slate-100 mb-2">
              Attribution Matrix
            </h2>
            <p className="text-sm text-slate-400">
              Analyze attribution matrix data
            </p>
          </Link>
        </div>
      </div>
    </div>
  );
}