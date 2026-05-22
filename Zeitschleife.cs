#region Copyright
///<remarks>
/// <Graz Lagrangian Particle Dispersion Model>
/// Copyright (C) [2019]  [Dietmar Oettl, Markus Kuntner]
/// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
/// the Free Software Foundation version 3 of the License
/// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
/// You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
///</remarks>
#endregion

using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading;

namespace GRAL_2001
{
    partial class Zeitschleife
    {
        private const float sqrt2F = 1.4142135623731F;
        private const float sqrtPiF = 1.77245385090552F;
        private const float aHurley = 0.1F;
        private const float bHurley = 0.6F;
        private const float a1Const = 0.05F;
        private const float a2Const = 1.7F;
        private const float a3Const = 1.1F;
        private const float Pi2F = 2F * MathF.PI;
        private const float RNG_Const = 2.328306435454494e-10F;
        private static float NextUniform01(ref uint m_z, ref uint m_w)
        {
            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
            uint u_rg = (m_z << 16) + m_w;
            return MathF.Max(1e-7f, (u_rg + 1) * RNG_Const);
        }

        private static float NextGaussianClamped(ref uint m_z, ref uint m_w)
        {
            float u1 = NextUniform01(ref m_z, ref m_w);
            float u2 = NextUniform01(ref m_z, ref m_w);

            float z = MathF.Sqrt(-2f * MathF.Log(u1)) * MathF.Sin(Pi2F * u2);
            return Math.Clamp(z, -2f, 2f);
        }

        // 20260510 Dense-stage entrainment random walk. Adds spatial dilution without removing pollutant mass.
        internal static void DenseHorzRandomStep(
            ref Program.DenseGasState state,
            float ux, float uy, float dt,
            ref uint m_z, ref uint m_w,
            out float denseRandX, out float denseRandY)
        {
            denseRandX = 0f;
            denseRandY = 0f;

            if (Program.GasBuoyancyMode != 1 || state.phase == Program.DensePhase.Neutral || !(dt > 0f))
            {
                return;
            }

            if (!Program.UseDenseHorzRandomWalk && !(Program.DenseEntrainmentFactor > 0f))
            {
                return;
            }

            float Hbulk = state.HbulkLocal;
            if (Hbulk <= 0f)
            {
                Hbulk = MathF.Max(Program.DenseGroundEps, state.cloudDepth);
            }
            Hbulk = MathF.Max(Program.DenseGroundEps, Hbulk);

            float KhDense = 0f;
            if (Program.UseDenseHorzRandomWalk && state.phase == Program.DensePhase.GroundSpread && state.Ug > 0f)
            {
                KhDense += MathF.Max(0f, Program.DenseHorzCoeff) * state.Ug * Hbulk;
            }

            float entrainment = MathF.Max(0f, Program.DenseEntrainmentFactor);
            if (entrainment > 0f)
            {
                float speedScale = MathF.Sqrt(MathF.Max(0f, ux * ux + uy * uy));
                speedScale = MathF.Max(speedScale, MathF.Abs(state.wp));
                speedScale = MathF.Max(speedScale, state.Ug);
                speedScale = MathF.Max(speedScale, 0.05f);

                KhDense += entrainment * speedScale * Hbulk;
            }

            if (!(KhDense > 0f))
            {
                return;
            }

            float sigmaDense = MathF.Sqrt(2f * KhDense * dt);
            float rxDense = NextGaussianClamped(ref m_z, ref m_w);
            float ryDense = NextGaussianClamped(ref m_z, ref m_w);

            float maxStep = MathF.Max(0f, Program.DenseHorzRandomMaxStepFrac) * MathF.Min(Program.GralDx, Program.DXK);
            if (maxStep > 0f)
            {
                denseRandX = Math.Clamp(sigmaDense * rxDense, -maxStep, maxStep);
                denseRandY = Math.Clamp(sigmaDense * ryDense, -maxStep, maxStep);
            }
            else
            {
                denseRandX = sigmaDense * rxDense;
                denseRandY = sigmaDense * ryDense;
            }
        }

        internal static void ApplyDenseLightBuoyancy(
            ref float ux, ref float uy, ref float uz,
            ref Program.DenseGasState state,
            float dt, float windSpeed,
            float heightAboveTerrain,
            float randSpread01)
        {
            if (Program.GasBuoyancyMode == 0) return;
            if (!(dt > 0f)) return;

            state.dbgDt = dt;
            state.dbgWind = windSpeed;

            float rhoAir = DenseRhoAir();
            bool useBulkState = Program.UseDenseGasMicroSpray;

            float rhoBulk;
            if (useBulkState)
            {
                if (!(state.rhoBulk > 0.01f))
                {
                    float rhoRel0 = Program.GasDensityRatio0;
                    if (!(rhoRel0 > 0f)) rhoRel0 = 1f;
                    rhoRel0 = MathF.Max(0.2f, MathF.Min(5f, rhoRel0));
                    state.rhoBulk = rhoAir * rhoRel0;
                }
                rhoBulk = MathF.Max(0.2f, state.rhoBulk);
            }
            else
            {
                float rhoRel0 = Program.GasDensityRatio0;
                if (!(rhoRel0 > 0f)) rhoRel0 = 1f;
                rhoRel0 = MathF.Max(0.2f, MathF.Min(5f, rhoRel0));
                rhoBulk = rhoAir * rhoRel0;
                state.rhoBulk = rhoBulk;
            }

            bool wantDense = (Program.GasBuoyancyMode == 1);
            float rhoRatio = rhoBulk / MathF.Max(1e-6f, rhoAir);
            float delta = MathF.Abs(rhoRatio - 1.0f);
            //bool denseAlive = state.BbulkLocal > Program.DenseEpsNeutral;//20260309
            float bFromLocal = state.BbulkLocal;
            float bFromRho = MathF.Max(0f, (rhoBulk - rhoAir) / MathF.Max(1e-6f, rhoAir));
            float bForGate = (bFromLocal > 0f) ? bFromLocal : bFromRho;
            bool denseAlive = bForGate > Program.DenseEpsNeutral;

            if (!wantDense)
            {
                const float g = 9.81f;
                float buoyAccelLight = g * (1.0f - rhoRatio);
                float uzMaxLight = 0.5f;
                uz += Math.Clamp(dt * buoyAccelLight * MathF.Max(0.0f, Program.DenseVertCoeff), -uzMaxLight, uzMaxLight);
                return;
            }

             if (!denseAlive)//20260309
            {
                state.phase = Program.DensePhase.Neutral;
                state.spreadActive = false;
                state.up = state.vp = state.wp = 0f;
                state.Ug = 0f;
                return;
            }

            if (heightAboveTerrain <= Program.DenseGroundEps)
            {
                if (state.phase != Program.DensePhase.GroundSpread)
                {
                    state.phase = Program.DensePhase.GroundSpread;
                    state.hitGroundOnce = true;
                    state.up = 0f;
                    state.vp = 0f;
                }
            }
            else
            {
                state.phase = Program.DensePhase.Descent;
            }

            // Keep existing descent logic for now.
            if (state.phase == Program.DensePhase.Descent)
            {
                const float g = 9.81f;
                float buoyAccel = g * (1.0f - rhoRatio);
                float vertScale = MathF.Max(0.0f, Program.DenseVertCoeff);
                float accel = buoyAccel * vertScale;

                if (useBulkState)
                {
                    float cDBase = 0.05f;
                    float cDWind = 0.10f * MathF.Min(1.0f, windSpeed / 5.0f);
                    float cDamp = MathF.Max(1e-4f, cDBase + cDWind);
                    float expd = MathF.Exp(-cDamp * dt);
                    state.wp = state.wp * expd + (accel / cDamp) * (1.0f - expd);

                    float neutralFactor = MathF.Max(0.0f, 1.0f - delta / 0.3f);
                    state.wp *= (1.0f - 0.5f * neutralFactor);

                    if (rhoRatio > 1.05f)
                    {
                        float hPull = 5.0f;
                        float wLow = MathF.Max(0.0f, 1.0f - heightAboveTerrain / hPull);
                        state.wp -= 0.1f * wLow * dt;
                    }

                    float wpMax = 0.5f;
                    state.wp = Math.Clamp(state.wp, -wpMax, wpMax);
                    uz += state.wp;

                    const float k0 = 0.0008f;
                    const float kU = 0.0004f;
                    const float kW = 0.0015f;
                    float k = MathF.Max(0.0f, Program.DenseMixFactor) * (k0 + kU * windSpeed + kW * MathF.Abs(state.wp));
                    float relax = 1.0f - MathF.Exp(-k * dt);

                    state.dbgWpAbs = MathF.Abs(state.wp);
                    state.dbgK = k;
                    state.dbgTau = 1f / MathF.Max(k, 1e-9f);
                    state.dbgRelaxCheck = 1f - MathF.Exp(-k0 * dt);
                    state.dbgRelax3s = 1f - MathF.Exp(-k * 3.0f);
                    state.dbgUgPhys = 0f;
                    state.dbgUgCapByStep = 0f;

                    state.lastRelax = relax;
                    state.rhoBulk = (1 - relax) * rhoBulk + relax * rhoAir;

                    if (MathF.Abs(state.rhoBulk - rhoAir) / MathF.Max(1e-6f, rhoAir) < Program.DenseEpsNeutral)
                    {
                        state.up = state.vp = state.wp = 0f;
                    }
                }
                else
                {
                    float uzMax = 0.5f;
                    uz += Math.Clamp(dt * accel, -uzMax, uzMax);
                }

                return;
            }

            // Ground spreading: use local column bulk from current horizontal cell.
            if (state.phase == Program.DensePhase.GroundSpread)
            {
                DenseEnsureSpreadDirection(ref state, randSpread01);
                //20260309
                // TODO: current GroundSpread intentionally ignores state.up/state.vp.
                // TODO: re-enable them only after Eq.(4a)-(4e) conservative momentum integration is implemented.
                float Hbulk = state.HbulkLocal;
                if (Hbulk <= 0f) Hbulk = MathF.Max(Program.DenseGroundEps, state.cloudDepth); // 工程回退：无局地列统计时用云团厚度

                float Bbulk = state.BbulkLocal;
                if (Bbulk <= 0f) Bbulk = MathF.Max(0f, (rhoBulk - rhoAir) / MathF.Max(1e-6f, rhoAir)); // 工程回退：用粒子体密度

                if (Bbulk < 0f) Bbulk = 0f;

                float UgPhys = Program.DenseSpreadFactor * MathF.Sqrt(MathF.Max(0f, 2f * 9.81f * Bbulk * Hbulk));
                float ugCapByStep = 0.5f * MathF.Min(Program.GralDx, Program.DXK) / MathF.Max(1e-3f, dt);
                state.Ug = MathF.Min(MathF.Min(UgPhys, Program.DenseUgMax), ugCapByStep);

                float Ugs = state.Ug * MathF.Cos(state.spreadAngle);
                float Vgs = state.Ug * MathF.Sin(state.spreadAngle);

                // only dense spreading velocity, no carry term
                ux += Ugs + state.up;
                uy += Vgs + state.vp;
                uz += state.wp;


                const float k0 = 0.0008f;
                const float kU = 0.0004f;
                const float kW = 0.0015f;
                float k = MathF.Max(0.0f, Program.DenseMixFactor) * (k0 + kU * windSpeed + kW * MathF.Abs(state.wp));
                float relax = 1.0f - MathF.Exp(-k * dt);

                state.dbgWpAbs = MathF.Abs(state.wp);
                state.dbgK = k;
                state.dbgTau = 1f / MathF.Max(k, 1e-9f);
                state.dbgRelaxCheck = 1f - MathF.Exp(-k0 * dt);
                state.dbgRelax3s = 1f - MathF.Exp(-k * 3.0f);

                state.lastRelax = relax;
                state.rhoBulk = (1 - relax) * rhoBulk + relax * rhoAir;

                if (MathF.Abs(state.rhoBulk - rhoAir) / MathF.Max(1e-6f, rhoAir) < Program.DenseEpsNeutral)
                {
                    state.phase = Program.DensePhase.Neutral;
                    state.spreadActive = false;
                    state.up = state.vp = state.wp = 0f;
                    state.Ug = 0f;
                }
            }
        }

        private static int DenseColumnLockIndex(int cellX, int cellY, int lockCount)
        {
            unchecked
            {
                int hash = cellX * 73856093 ^ cellY * 19349663;
                if (hash < 0) hash = -hash;
                return hash % Math.Max(1, lockCount);
            }
        }

        private static void DenseColumnAccumulate(int cellX, int cellY, float heightAboveTerrain, float rhoParticle, float rhoAir)
        {
            if (Program.DenseColumnCount == null || Program.DenseColumnCount.Length == 0) return;
            if (cellX < 0 || cellX >= Program.DenseColumnCount.Length) return;
            if (cellY < 0 || cellY >= Program.DenseColumnCount[cellX].Length) return;
            if (Program.DenseColumnLocks == null || Program.DenseColumnLocks.Length == 0) return;

            int lockId = DenseColumnLockIndex(cellX, cellY, Program.DenseColumnLocks.Length);
            lock (Program.DenseColumnLocks[lockId])
            {
                Program.DenseColumnSumH[cellX][cellY] += heightAboveTerrain;
                Program.DenseColumnSumRho[cellX][cellY] += rhoParticle;
                int count = ++Program.DenseColumnCount[cellX][cellY];

                float hbulk = (float)(Program.DenseColumnSumH[cellX][cellY] / Math.Max(1, count));
                float rhobulk = (float)(Program.DenseColumnSumRho[cellX][cellY] / Math.Max(1, count));
                float bbulk = MathF.Max(0f, (rhobulk - rhoAir) / MathF.Max(1e-6f, rhoAir));

                Program.DenseColumnHbulk[cellX][cellY] = hbulk;
                Program.DenseColumnRhoBulk[cellX][cellY] = rhobulk;
                Program.DenseColumnBbulk[cellX][cellY] = bbulk;
            }
        }

        private static void DenseColumnSample(int cellX, int cellY, float rhoAir, out float hbulk, out float rhobulk, out float bbulk, out int count)
        {
            hbulk = 0f;
            rhobulk = rhoAir;
            bbulk = 0f;
            count = 0;

            if (Program.DenseColumnCount == null || Program.DenseColumnCount.Length == 0) return;
            if (cellX < 0 || cellX >= Program.DenseColumnCount.Length) return;
            if (cellY < 0 || cellY >= Program.DenseColumnCount[cellX].Length) return;

            if (!Program.DenseColumnSmooth3x3)
            {
                count = Program.DenseColumnCount[cellX][cellY];
                if (count > 0)
                {
                    hbulk = Program.DenseColumnHbulk[cellX][cellY];
                    rhobulk = Program.DenseColumnRhoBulk[cellX][cellY];
                    bbulk = Program.DenseColumnBbulk[cellX][cellY];
                }
                return;
            }

            double sumH = 0.0;
            double sumRho = 0.0;
            int sumCount = 0;

            int iMin = Math.Max(1, cellX - 1);
            int iMax = Math.Min(Program.DenseColumnCount.Length - 2, cellX + 1);
            int jMin = Math.Max(1, cellY - 1);
            int jMax = Math.Min(Program.DenseColumnCount[cellX].Length - 2, cellY + 1);

            for (int i = iMin; i <= iMax; i++)
            {
                for (int j = jMin; j <= jMax; j++)
                {
                    int c = Program.DenseColumnCount[i][j];
                    if (c <= 0) continue;
                    sumH += Program.DenseColumnHbulk[i][j] * c;
                    sumRho += Program.DenseColumnRhoBulk[i][j] * c;
                    sumCount += c;
                }
            }

            count = sumCount;
            if (sumCount <= 0) return;

            hbulk = (float)(sumH / sumCount);
            rhobulk = (float)(sumRho / sumCount);
            bbulk = MathF.Max(0f, (rhobulk - rhoAir) / MathF.Max(1e-6f, rhoAir));
        }

        private static float DenseRhoAir()
        {
            return 1.2041f * (273.15f / MathF.Max(200f, Program.AmbientTempK));
        }

        // Eq(6b): wr = wi * 0.5 * (tanh(-9 * rhoRatio + 11.3) + 1)
        private static float DenseReflectFactor(float rhoRatio)
        {
            float f = 0.5f * (MathF.Tanh(-9.0f * rhoRatio + 11.3f) + 1.0f);
            return Math.Clamp(f, 0f, 1f);
        }
        
        //新增统一反射 helper 20260312
        private static void ApplyDenseBottomReflection(ref float velzold, ref Program.DenseGasState st)
        {
            // Neutral particles: use default elastic reflection
            if (st.phase == Program.DensePhase.Neutral)
            {
                velzold = -velzold;
                st.dbgWr = MathF.Abs(velzold);
                return;
            }

            float rhoAir = DenseRhoAir();
            float rhoBulkForRefl = (st.RhoBulkLocal > 0f) ? st.RhoBulkLocal : st.rhoBulk;
            float rhoRatio = rhoBulkForRefl / MathF.Max(1e-6f, rhoAir);

            // Incident speed toward ground; fallback keeps robustness
            float wi = MathF.Max(0f, -velzold);
            if (wi <= 0f) wi = MathF.Abs(velzold);

            float f = DenseReflectFactor(rhoRatio); // Eq(6b)
            float wr = wi * f;                      // Eq(6a)

            st.hitGroundOnce = true;
            if (st.phase != Program.DensePhase.Neutral)
            {
                st.phase = Program.DensePhase.GroundSpread;
            }
            st.dbgWr = wr;

            // reflected upward speed
            velzold = wr;
        }

        private static void DenseEnsureSpreadDirection(ref Program.DenseGasState state, float rand01)
        {
            if (!state.spreadActive)
            {
                state.spreadAngle = MathF.PI * 2f * Math.Clamp(rand01, 0f, 0.999999f);
                state.spreadActive = true;
            }
        }


        /// <summary>
        ///Time loop for all released particles - the particles are tracked until they leave the domain area(steady state mode)
        ///or until the dispersion time has expired (transient mode) 
        /// </summary>
        /// <param name="nteil">Particle number</param>
        /// /// <remarks>
        /// This class is the particle driver for particles, released by a source in steady state and transient GRAL mode
        /// </remarks>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void Calculate(int nteil)
        {
            bool endOfStep = false;//20260113
            //random number generator seeds
            int rnd = (Environment.TickCount + nteil) & Int32.MaxValue;
            uint m_w = (uint)(rnd + 521288629);
            uint m_z = (uint)(rnd + 2232121);
            if (Program.UseFixedRndSeedVal)
            {
                m_w = Program.RnGSeed.Seed1 + (uint)nteil;
                m_z = Program.RnGSeed.Seed2 + (uint)nteil * 2;
            }

            float zahl1 = 0;
            uint u_rg = 0;
            float u1_rg = 0;

            //local particle coordinates and pollution "mass"
            double xcoord_nteil = Program.Xcoord[nteil];
            double ycoord_nteil = Program.YCoord[nteil];
            float zcoord_nteil = Program.ZCoord[nteil];
            double masse = Program.ParticleMass[nteil];
            if (nteil <= 3 && Program.IWET == 1)//20260301
            {
                System.Console.WriteLine($"[DBG] init n={nteil} x={xcoord_nteil} y={ycoord_nteil} z={zcoord_nteil}");
            }

            //get index of internal source group number
            int SG_nteil = Program.ParticleSG[nteil];
            for (int i = 0; i < Program.SourceGroups.Count; i++)
            {
                if (SG_nteil == Program.SourceGroups[i])
                {
                    SG_nteil = i;
                    break;
                }
            }

            // deposition parameters
            int Deposition_type = Program.ParticleMode[nteil]; // 0 = no deposition, 1 = depo + conc, 2 = only deposition
            float vsed = Program.ParticleVsed[nteil];
            float vdep = Program.ParticleVdep[nteil];
            double area_rez_fac = Program.TAUS * Program.GridVolume * 24 / 1000; // conversion factor for particle mass from 碌g/m鲁/s to mg per day
            area_rez_fac = area_rez_fac / (Program.GralDx * Program.GralDy);  // conversion to mg/m虏

            double xtun = 0;
            double ytun = 0;

            // variables for the Hurley plume rise algorithm
            float MeanHurley = 0.1F; // start value
            float GHurley = 0;
            float FHurley = 0;
            float MHurley = 0;
            float RHurley = 0;
            float wpHurley = 0;
            float upHurley = 0;
            float DeltaZHurley = 0;
            float SigmaUpHurley = 0;
            float minVDI3782 = 0; //Minimum plume rise height for cold stacks depending on VDI 3782-3 Equation 18
            //float DeltaZHurleySumme = 0;

            //Grid variables
            int GrammCellX = 1, GrammCellY = 1; // default value for flat terrain
            int FFCellX, FFCellY;               // Flow field cells
            float FFGridX = Program.DXK;
            float FFGridY = Program.DYK;
            float FFGridXRez = 1 / Program.DXK;
            float FFGridYRez = 1 / Program.DYK;
            float GrammGridXRez = 1 / Program.DDX[1];
            float GrammGridYRez = 1 / Program.DDY[1];
            float ConcGridXRez = 1 / Program.GralDx;
            float ConcGridYRez = 1 / Program.GralDy;
            float ConcGridZRez = 2 / Program.GralDz;
            float ConcGridXHalf = Program.GralDx * 0.5F;
            float ConcGridYHalf = Program.GralDy * 0.5F;

            int IFQ = Program.TS_Count;
            int Kenn_NTeil = Program.ParticleSource[nteil];       // number of source
            if (Kenn_NTeil == 0)
            {
                return;
            }

            int SourceType = Program.SourceType[nteil]; // sourcetype
            float blh = Program.BdLayHeight;
            int IKOOAGRAL = Program.IKOOAGRAL;
            int JKOOAGRAL = Program.JKOOAGRAL;

            Span<int> kko = new int[Program.NS];
            Span<double> ReceptorConcentration = new double[Program.ReceptorNumber + 1];

            int reflexion_flag = Consts.ParticleNotReflected;
            int ISTATISTIK = Program.IStatistics;
            int ISTATIONAER = Program.ISTATIONAER;
            int topo = Program.Topo;

            int TransientGridX = 1;
            int TransientGridY = 1;
            int TransientGridZ = 1;
            Boolean emission_timeseries = Program.EmissionTimeseriesExist;

            double Area_cart = Program.DXK * Program.DYK;

            int reflexion_number = 0; 		  // counter to limit max. number of reflexions
            int timestep_number = 0;         // counter for the time-steps LOG-output

            double decay_rate = Program.DecayRate[Program.ParticleSG[nteil]]; // decay rate for real source group number of this source

            //for transient simulations dispersion time needs to be equally distributed from zero to TAUS
            float tgesamt = Program.DispTimeSum;
            if (ISTATIONAER == Consts.TransientMode)
            {
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                uint u9 = (m_z << 16) + m_w;
                float zuff1 = (float)((u9 + 1.0) * 2.328306435454494e-10);
                tgesamt = zuff1 * tgesamt;

                //modulate emission
                if (emission_timeseries == true && (!Program.ContinuousTraj || (!Program.TransientReleaseEnabled && Program.IWET == 1))) //20260521
                //{
                    //if (Program.IWET <= Program.EmFacTimeSeries.GetUpperBound(0) && SG_nteil <= Program.EmFacTimeSeries.GetUpperBound(1))
                    //{
                    //    masse *= Program.EmFacTimeSeries[Program.IWET - 1, SG_nteil];
                    //}
                //if (emission_timeseries)
                {
                    float fac = 0f;
                    int row = Program.IWET - 1;

                    if (row >= 0 &&
                        row <= Program.EmFacTimeSeries.GetUpperBound(0) &&
                        SG_nteil <= Program.EmFacTimeSeries.GetUpperBound(1))
                    {
                        fac = Program.EmFacTimeSeries[row, SG_nteil];
                    }

                    masse *= fac; 

                    if (masse <= 0)
                    {
                        Program.Remove_MassZero++;//20260125
                        goto REMOVE_PARTICLE;
                    }
                }
                //try
                if (Deposition_type < Consts.DepoOnly) // compute concentrations for this particle
                {
                    lock (Program.EmissionPerSG)
                    {
                        Program.EmissionPerSG[SG_nteil] += masse;
                    }
                }
                //catch{}
            }
          

            /*
             *   INITIIALIZING PARTICLE PROPERTIES -> DEPEND ON SOURCE CATEGORY
             */

            float AHint = 0;
            float AHintold = 0;
            float varw = 0;
            double xturb = 0;  // entrainment turbulence for tunnel jet streams in x-direction
            double yturb = 0;  // entrainment turbulence for tunnel jet streams in y-direction

            //advection time step
            float idt = 0.1F;
            float auszeit = 0;

            //flag indicating if particle belongs to a tunnel portal 1 = not a portal
            int tunfak = Consts.ParticleIsNotAPortal;

            float DSIGUDX = 0, DSIGUDY = 0, DSIGVDX = 0, DSIGVDY = 0;
            float DUDX = 0, DVDX = 0, DUDY = 0, DVDY = 0;

            //relative coordinates of particles inside model domain (used to calc. cell-indices of particles moving in the microscale flow field)
            double xsi = xcoord_nteil - IKOOAGRAL;
            double eta = ycoord_nteil - JKOOAGRAL;
            double xsi1 = xcoord_nteil - Program.GrammWest;
            double eta1 = ycoord_nteil - Program.GrammSouth;
            if ((eta <= Program.EtaMinGral) || (xsi <= Program.XsiMinGral) || (eta >= Program.EtaMaxGral) || (xsi >= Program.XsiMaxGral))
            {
                if (Program.Remove_OutDomain < 10)
                {
                    string outreason =
                        (xsi <= Program.XsiMinGral) ? "xsi<=min" :
                        (xsi >= Program.XsiMaxGral) ? "xsi>=max" :
                        (eta <= Program.EtaMinGral) ? "eta<=min" :
                        "eta>=max";

                    System.Console.WriteLine(
                        $"[DBG] OUT[{outreason}] n={nteil} IWET={Program.IWET} " +
                        $"x={xcoord_nteil} y={ycoord_nteil} z={zcoord_nteil} " +
                        $"PartH={(zcoord_nteil - AHint)} AHint={AHint} " +
                        $"xsi={xsi} eta={eta} XsiMin={Program.XsiMinGral} XsiMax={Program.XsiMaxGral} " +
                        $"EtaMin={Program.EtaMinGral} EtaMax={Program.EtaMaxGral}"
                    );
                }                
                Program.Remove_OutDomain++;//20260114
                goto REMOVE_PARTICLE;
            }

            int Max_Loops = (int)(3E6 + Math.Min(9.7E7,
                Math.Max(Math.Abs(Program.EtaMaxGral - Program.EtaMinGral), Math.Abs(Program.XsiMaxGral - Program.XsiMinGral)) * 1000)); // max. Loops 1E6 (max. nr. of reflexions) + 2E6 Min + max(x,y)/0.001 
            int Max_Reflections = Math.Min(1000000, Program.NII * Program.NJJ * 10); // max. 10 reflections per cell
            int ConcCellPrevX = -1; int ConcCellPrevY = -1; int ConcCellPrevZ = 0; float ConcCellTimeMax = Math.Max(100, Program.GralDx / 0.04F); float ConcCellTime = 0; float distanceParticle = 50;
            
            //interpolated orography
            float PartHeightAboveTerrain = 0;
            float PartHeightAboveBuilding = 0;

            FFCellX = (int)(xsi * FFGridXRez) + 1;
            FFCellY = (int)(eta * FFGridYRez) + 1;
            if ((FFCellX < 1) || (FFCellX > Program.NII) || (FFCellY < 1) || (FFCellY > Program.NJJ))
            {
                goto REMOVE_PARTICLE;
            }

            if (topo == Consts.TerrainAvailable)
            {
                //with terrain
                GrammCellX = Math.Clamp((int)(xsi1 * GrammGridXRez) + 1, 1, Program.NX);
                GrammCellY = Math.Clamp((int)(eta1 * GrammGridYRez) + 1, 1, Program.NY);

                //interpolated orography
                AHint = Program.AHK[FFCellX][FFCellY];
                PartHeightAboveTerrain = zcoord_nteil - AHint;
                PartHeightAboveBuilding = PartHeightAboveTerrain;
            }
            else
            {
                AHint = 0;  
                //flat terrain
                PartHeightAboveTerrain = zcoord_nteil;
                PartHeightAboveBuilding = PartHeightAboveTerrain;
            }

            //wind-field interpolation
            float UXint = 0, UYint = 0, UZint = 0;
            int IndexK = 1;
            (UXint, UYint, UZint, IndexK) = IntWindCalculate(FFCellX, FFCellY, AHint, xcoord_nteil, ycoord_nteil, zcoord_nteil);
            float windge = MathF.Max(0.1F, MathF.Sqrt(Program.Pow2(UXint) + Program.Pow2(UYint)));
            if (float.IsNaN(UXint) || float.IsNaN(UYint) || float.IsNaN(UZint) ||
            float.IsInfinity(UXint) || float.IsInfinity(UYint) || float.IsInfinity(UZint))//20260114
            {
                UXint = 0f; UYint = 0f; UZint = 0f;
            }

            //variables tunpa needs to be defined even if not used
            float tunpa = 0;

            float ObL = Program.Ob[GrammCellX][GrammCellY];
            float roughZ0 = Program.Z0Gramm[GrammCellX][GrammCellY];
            if (Program.AdaptiveRoughnessMax > 0)
            {
                ObL = Program.OLGral[FFCellX][FFCellY];
                roughZ0 = Program.Z0Gral[FFCellX][FFCellY];
            }
            float deltaZ = roughZ0;

            //vertical interpolation of horizontal standard deviations of wind component fluctuations between observations
            (float U0int, float V0int) = IntStandCalculate(nteil, roughZ0, PartHeightAboveBuilding, windge, SigmaUpHurley);

            //remove particles above boundary-layer
            if ((PartHeightAboveBuilding > blh) && (Program.ISTATIONAER != Consts.TransientMode) && (ObL < 0)) //26042020 (Ku): removed for transient mode -> particles shoulb be tracked above blh
            {
                goto REMOVE_PARTICLE;
            }

            //initial turbulent velocities
            float velxold = 0;
            float velyold = 0;
            float velzold = 0;

            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
            u_rg = (m_z << 16) + m_w;
            u1_rg = (u_rg + 1) * RNG_Const;
            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
            u_rg = (m_z << 16) + m_w;
            zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

            velxold = U0int * zahl1;

            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
            u_rg = (m_z << 16) + m_w;
            u1_rg = (u_rg + 1) * RNG_Const;
            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
            u_rg = (m_z << 16) + m_w;
            zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

            velyold = V0int * zahl1;

            //dispersion time of a particle
            auszeit = 0;

            //initial properties for particles stemming from point sources
            bool useExplicitExitVelocityDirection = false; //20260521 1852 per-particle switch: explicit vector initializes velocity once and disables Hurley rise
            if (SourceType == Consts.SourceTypePoint)
            {
                float ExitVelocity = Program.PS_V[Kenn_NTeil];
                int velTimeSeriesIndex = Program.PS_TimeSeriesVelocity[Kenn_NTeil];
                if (velTimeSeriesIndex != -1) // Time series for exit velocity
                {
                    if ((Program.IWET - 1) < Program.PS_TimeSerVelValues[velTimeSeriesIndex].Value.Length)
                    {
                        ExitVelocity = Program.PS_TimeSerVelValues[velTimeSeriesIndex].Value[Program.IWET - 1];
                    }
                }
                ExitVelocity = MathF.Max(0, ExitVelocity);

                // direction of jet velocity 20260408; gated by cfg switch 20260521 1852
                float dirx = Program.PS_VDirX[Kenn_NTeil];
                float diry = Program.PS_VDirY[Kenn_NTeil];
                float dirz = Program.PS_VDirZ[Kenn_NTeil];

                float dnorm = MathF.Sqrt(dirx * dirx + diry * diry + dirz * dirz);
                useExplicitExitVelocityDirection = Program.UseExplicitExitVelocityDirection && dnorm > 1e-6f; //20260521 1852 Trial 8 opt-in only
                if (useExplicitExitVelocityDirection)
                {
                    dirx /= dnorm;
                    diry /= dnorm;
                    dirz /= dnorm;
                }
                else
                {
                    dirx = 0f; diry = 0f; dirz = 1f;
                }

                float jetVx = ExitVelocity * dirx;
                float jetVy = ExitVelocity * diry;
                float jetVz = ExitVelocity * dirz;

                if (useExplicitExitVelocityDirection) //20260521 1852 initialize explicit exit velocity once at particle creation
                {
                    velxold += jetVx;
                    velyold += jetVy;
                    velzold += jetVz;
                }
                else if (jetVz < 0f) //20260521 1852 keep legacy downward correction when explicit vector mode is disabled
                {
                    velzold += jetVz;
                }

                float ExitVelocityUp = useExplicitExitVelocityDirection ? 0f : MathF.Max(0f, jetVz); //20260521 1852 explicit vector mode owns momentum, so Hurley gets no upward velocity
                
                float ExitTemperature = Program.PS_T[Kenn_NTeil];
                int tempTimeSeriesIndex = Program.PS_TimeSeriesTemperature[Kenn_NTeil];
                if (tempTimeSeriesIndex != -1) // Time series for exit temperature
                {
                    if ((Program.IWET - 1) < Program.PS_TimeSerTempValues[tempTimeSeriesIndex].Value.Length)
                    {
                        ExitTemperature = Program.PS_TimeSerTempValues[tempTimeSeriesIndex].Value[Program.IWET - 1] + 273;
                    }
                }
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                u1_rg = (u_rg + 1) * RNG_Const;
                
                //fluctuation of exit / ambient temperature
                ExitTemperature = MathF.Max(273, ExitTemperature * (1.05F - u1_rg * 0.2F));

                //plume-rise velocity
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                u1_rg = (u_rg + 1) * RNG_Const;
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);
                                                
                FHurley = (9.81F * ExitVelocityUp * Program.Pow2(Program.PS_D[Kenn_NTeil] * 0.5F) *
                                 (ExitTemperature - 273F) / ExitTemperature);
                GHurley = (273 / ExitTemperature * ExitVelocityUp * Program.Pow2(Program.PS_D[Kenn_NTeil] * 0.5F));
                MHurley = GHurley * ExitVelocityUp;
                RHurley = MathF.Sqrt(ExitVelocityUp / MathF.Sqrt(Program.Pow2(windge) + Program.Pow2(ExitVelocity)));
                wpHurley = ExitVelocityUp;
                upHurley = MathF.Sqrt(Program.Pow2(windge) + Program.Pow2(wpHurley));
                MeanHurley = 0;
                DeltaZHurley = 0;
                //minVDI3782 = 3 * ExitVelocityup * Program.PS_D[Kenn_NTeil] / windge; //Minimum plume rise height for cold stacks depending on VDI 3782-3
                minVDI3782 = 3 * ExitVelocityUp * Program.PS_D[Kenn_NTeil] / MathF.Max(windge, 0.1F);

            }

            //initial properties for particles stemming from tunnel portals
            double vtunx = 0; double vtuny = 0; float tunbe = 0; double tuncos = 0; double tunsin = 0; double tunqu = 0;
            double tunx = 0; double tuny = 0; float tunbreitalt = 0; float tunbreit = 0; float yteilchen = 0;
            double yneu = 0; double xneuy = 0; double yneuy = 0; float zturb = 0; double yold = 0;
            float coruri = 0; float corqri = 0; float coswind = 0; float sinwind = 0; float UXneu = 0;
            float UYneu = 0; float xrikor = 0; float yrikor = 0; float buoy = 0; float B = 0;
            float HydrD = 0; float TLT0 = 0; /*float constA = 0.35F;*/ float faktzeit = 1; float geschwalt = 0; float geschwneu = 0;

            float TS_ExitTemperature = 0;
            float TS_ExitVelocity = 0;
            if (SourceType == Consts.SourceTypePortal) // Tunnel portals
            {
                TS_ExitVelocity = Program.TS_V[Kenn_NTeil];
                int velTimeSeriesIndex = Program.TS_TimeSeriesVelocity[Kenn_NTeil];
                if (velTimeSeriesIndex != -1) // Time series for exit velocity
                {
                    if ((Program.IWET - 1) < Program.TS_TimeSerVelValues[velTimeSeriesIndex].Value.Length)
                    {
                        TS_ExitVelocity = Program.TS_TimeSerVelValues[velTimeSeriesIndex].Value[Program.IWET - 1];
                    }
                }

                TS_ExitTemperature = Program.TS_T[Kenn_NTeil];
                int tempTimeSeriesIndex = Program.TS_TimeSeriesTemperature[Kenn_NTeil];
                if (tempTimeSeriesIndex != -1) // Time series for exit temperature
                {
                    if ((Program.IWET - 1) < Program.TS_TimeSerTempValues[tempTimeSeriesIndex].Value.Length)
                    {
                        TS_ExitTemperature = Program.TS_TimeSerTempValues[tempTimeSeriesIndex].Value[Program.IWET - 1];
                    }
                }

                vtunx = TS_ExitVelocity * (Program.TS_Y2[Kenn_NTeil] - Program.TS_Y1[Kenn_NTeil]) / Program.TS_Width[Kenn_NTeil];
                vtuny = TS_ExitVelocity * (Program.TS_X1[Kenn_NTeil] - Program.TS_X2[Kenn_NTeil]) / Program.TS_Width[Kenn_NTeil];
                tunbe = (float)Math.Sqrt(Program.Pow2(vtunx) + Program.Pow2(vtuny));
                tuncos = vtunx / (tunbe + 0.01F);
                tunsin = vtuny / (tunbe + 0.01F);
                tunqu = 0;
                tunpa = tunbe;
                tunfak = Consts.ParticleIsAPortal;
                tunx = tunpa * tuncos - tunqu * tunsin;
                tuny = tunpa * tunsin + tunqu * tuncos;

                tunbreitalt = Program.FloatMax(Program.TS_Width[Kenn_NTeil], 0.00001F);
                tunbreit = tunbreitalt;
                xtun = 0.5F * (Program.TS_X1[Kenn_NTeil] + Program.TS_X2[Kenn_NTeil]);
                ytun = 0.5F * (Program.TS_Y1[Kenn_NTeil] + Program.TS_Y2[Kenn_NTeil]);

                //relative particle position to the jet stream centre-line
                yteilchen = (float)(-xcoord_nteil * tunsin + ycoord_nteil * tuncos + xtun * tunsin - ytun * tuncos);
                yold = yteilchen;
                yneu = yold;
                xneuy = xtun - yteilchen * tunsin;
                yneuy = ytun + yteilchen * tuncos;
                zturb = 0;

                //influence of changing wind directions on the position of the jet-stream centre-line
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                u1_rg = (u_rg + 1) * RNG_Const;
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

                coruri = windge + U0int * (float) zahl1;

                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                u1_rg = (u_rg + 1) * RNG_Const;
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

                corqri = V0int * zahl1;

                coswind = UXint / windge;
                sinwind = UYint / windge;
                UXneu = coruri * coswind - corqri * sinwind;
                UYneu = coruri * sinwind + corqri * coswind;
                if (UXint == 0)
                {
                    UXint = 0.001F;
                }

                if (UYint == 0)
                {
                    UYint = 0.001F;
                }

                xrikor = UXneu / UXint;
                yrikor = UYneu / UYint;

                //bouyancy of the tunnel jet
                buoy = 0;
                B = TS_ExitTemperature / 283 * 9.81F;
                HydrD = 4 * Program.TS_Height[Kenn_NTeil] * Program.TS_Width[Kenn_NTeil] / (2 * (Program.TS_Height[Kenn_NTeil] + Program.TS_Width[Kenn_NTeil]));
                if (B >= 0)
                {
                    TLT0 = MathF.Sqrt(HydrD / (1.8F * (B + 0.00001F)));
                }
                else
                {
                    B = 0;
                }
                //double dummy = corqri + coruri + coswind + sinwind + UXneu + UYneu + xrikor + yrikor + B + HydrD + TLT0;
                AHintold = AHint;
            }

            int depo_reflection_counter = -5; // counter to ensure, a particle moves top to down
            bool TerrainStepAllowed = true;   //19.05.25 Ku: Flag, that one step is allowed
            int depo_events = 0; // deposition counter
            int ground_hit_events = 0; // ground-hit counter
        /*
        *      LOOP OVER THE TIME-STEPS
        */
            float ObLength = 0;
            float Ustern = 0;
 
        MOVE_FORWARD:
            while (timestep_number <= Max_Loops)
            {
                ++timestep_number;
                double xcoord_nteil_Prev = xcoord_nteil;
                double ycoord_nteil_Prev = ycoord_nteil;
                float zcoord_nteil_Prev = zcoord_nteil;
                int FFCellXPrev = FFCellX;
                int FFCellYPrev = FFCellY;

                //if particle is within the user-defined tunnel-entrance zone it is removed (sucked-into the tunnel)
                if (Program.TunnelEntr == true)
                {
                    if ((Program.TUN_ENTR[FFCellX][FFCellY] == 1) && (PartHeightAboveTerrain <= 5))
                    {
                        Program.Remove_Tunnel++; // 涓存椂璁℃暟锛歍unnelEntr-20260114
                        goto REMOVE_PARTICLE;
                    }
                }

                //interpolate wind field
                (UXint, UYint, UZint, IndexK) = IntWindCalculate(FFCellX, FFCellY, AHint, xcoord_nteil, ycoord_nteil, zcoord_nteil);
                windge = MathF.Max(0.01F, MathF.Sqrt(Program.Pow2(UXint) + Program.Pow2(UYint)));
                if (float.IsNaN(windge) || float.IsInfinity(windge)) windge = 0.01f;//20260114
                if (float.IsNaN(UXint) || float.IsNaN(UYint) || float.IsNaN(UZint) ||
                float.IsInfinity(UXint) || float.IsInfinity(UYint) || float.IsInfinity(UZint))//20260114
                {
                    UXint = 0f; UYint = 0f; UZint = 0f;
                }

                //particles below the surface or buildings are removed
                if (PartHeightAboveTerrain < 0)
                {
                    ground_hit_events++; // ground-hit counter
                    if (Program.ContinuousTraj) // keep particle in continuous trajectory mode
                    {
                        Program.Keep_Ground++;
                        if (Program.Keep_Ground < 5)
                        {
                            System.Console.WriteLine($"[DBG] KEEP_GROUND z={zcoord_nteil} AHint={AHint}");
                        }

                        PartHeightAboveTerrain = 0;
                        zcoord_nteil = AHint + 0.01f;

                        if (Program.GasBuoyancyMode == 1 && Program.DenseStates != null &&
                            nteil > 0 && nteil < Program.DenseStates.Length)
                        {
                            ref Program.DenseGasState st = ref Program.DenseStates[nteil];
                            st.hitGroundOnce = true;
                            if (st.phase != Program.DensePhase.Neutral)
                            {
                                st.phase = Program.DensePhase.GroundSpread;
                            }
                        }

                        reflexion_flag = Consts.ParticleReflected;
                    }
                    else
                    {
                        Program.Remove_Ground++;
                        goto REMOVE_PARTICLE;
                    }
                }

                /*
                 *   VERTICAL DIFFUSION ACCORING TO FRANZESE, 1999
                 */

                //dissipation            
                if (Program.AdaptiveRoughnessMax > 0)
                {
                    ObLength = Program.OLGral[FFCellX][FFCellY];
                    Ustern = Program.USternGral[FFCellX][FFCellY];
                    roughZ0 = Program.Z0Gral[FFCellX][FFCellY];
                }
                else
                {
                    ObLength = Program.Ob[GrammCellX][GrammCellY];
                    Ustern = Program.Ustern[GrammCellX][GrammCellY];
                    roughZ0 = Program.Z0Gramm[GrammCellX][GrammCellY];
                }
                deltaZ = roughZ0;
                if (deltaZ < PartHeightAboveBuilding)
                {
                    deltaZ = PartHeightAboveBuilding;
                }

                float eps = Program.Pow3(Ustern) / deltaZ / 0.4F *
                                    MathF.Pow(1F + 0.5F * MathF.Pow(PartHeightAboveBuilding / Math.Abs(ObLength), 0.8F), 1.8F);
                if (float.IsNaN(Ustern) || float.IsInfinity(Ustern)) Ustern = 0.15f;
                if (float.IsNaN(ObLength) || float.IsInfinity(ObLength)) ObLength = 1000f;//20260114
                if (float.IsNaN(roughZ0) || float.IsInfinity(roughZ0)) roughZ0 = 0.01f;//20260114
                if (float.IsNaN(PartHeightAboveBuilding) || float.IsInfinity(PartHeightAboveBuilding)) PartHeightAboveBuilding = 0.01f;//20260114
                if (float.IsNaN(eps) || float.IsInfinity(eps) || eps < 1e-6f)//20260114
                {
                    eps = 1e-6f;
                }

                depo_reflection_counter++;
                float W0int = 0; float varw2 = 0; float skew = 0; float curt = 0; float dskew = 0;
                float wstern = Ustern * 2.9F;
                float vert_blh = PartHeightAboveBuilding / blh;
                float hterm = (1 - vert_blh);
                /*float varwu = 0;*/ float dvarw = 0; float dcurt = 0;

                if (ObLength >= 0)
                {
                    if (ISTATISTIK != Consts.MeteoSonic)
                    {
                        W0int = Ustern * Program.StdDeviationW * 1.25F;
                    }
                    else
                    {
                        W0int = Program.W0int;
                    }

                    varw = Program.Pow2(W0int);
                    varw2 = Program.Pow2(varw);
                    curt = 3 * varw2;
                }
                else
                {
                    if (ISTATISTIK == Consts.MeteoSonic)
                    {
                        varw = Program.Pow2(Program.W0int);
                    }
                    else
                    {
                        varw = Program.Pow2(Ustern) * Program.Pow2(1.15F + 0.1F * MathF.Pow(blh / (-ObLength), 0.67F)) * Program.Pow2(Program.StdDeviationW);
                    }
                    varw2 = Program.Pow2(varw);

                    if ((PartHeightAboveBuilding < Program.BdLayH10) || (PartHeightAboveBuilding >= blh))
                    {
                        skew = 0;
                        dskew = 0;
                    }
                    else
                    {
                        skew = a3Const * PartHeightAboveBuilding * Program.Pow2(hterm) * Program.Pow3(wstern) / blh;
                        dskew = a3Const * Program.Pow3(wstern) / blh * Program.Pow2(hterm) - 2 * a3Const * Program.Pow3(wstern) * PartHeightAboveBuilding * hterm / Program.Pow2(blh);
                    }

                    if (PartHeightAboveBuilding < Program.BdLayH10)
                    {
                        curt = 3 * varw2;
                    }
                    else
                    {
                        curt = 3.5F * varw2;
                    }
                }

                //time-step calculation
                if (tunfak == Consts.ParticleIsNotAPortal)
                {
                    idt = 0.1F * varw / eps;
                    //in case plume rise the time step shall not exceed 0.2s in the first 20 seconds
                    if (SourceType == Consts.SourceTypePoint && auszeit < 20)
                    {
                        idt = MathF.Min(0.2F, idt);
                    }

                    float idtt = 0.5F * Math.Min(Program.GralDx, FFGridX) / (windge + MathF.Sqrt(Program.Pow2(velxold) + Program.Pow2(velyold)));
                    if (idtt < idt)
                    {
                        idt = idtt;
                    }

                    if (float.IsNaN(idt) || float.IsInfinity(idt)) idt = 0.1f;//20260202
                    idt = Math.Clamp(idt, 0.01F, 0.5F);
                }

                float dummyterm = curt - Program.Pow2(skew) / varw - varw2;
                if (dummyterm == 0)
                {
                    dummyterm = 0.000001F;
                }

                float alpha = (0.3333F * dcurt - skew / (2 * varw) * (dskew - Program.C0z * eps) - varw * dvarw) / dummyterm;
                //alpha needs to be bounded to keep algorithm stable
                Math.Clamp(alpha, -0.1F, 0.1F);

                float beta = (dskew - 2 * skew * alpha - Program.C0z * eps) / (2 * varw);
                float gamma = dvarw - varw * alpha;

                //particle acceleration
                float acc = alpha * Program.Pow2(velzold) + beta * velzold + gamma;

                //new vertical wind speed
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                u1_rg = (u_rg + 1) * RNG_Const;
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);
                Math.Clamp(zahl1, -2, 2);
                float velz = acc * idt + MathF.Sqrt(Program.C0z * eps * idt) * zahl1 + velzold;
                
                //******************************************************************************************************************** OETTL, 31 AUG 2016
                //in adjecent cells to vertical solid walls, turbulent velocities are only allowed in the direction away from the wall
                if (velz == 0)
                {
                    velz = 0.01F;
                }

                short _kkart = Program.KKART[FFCellX][FFCellY];
                //solid wall west of the particle
                if ((Program.KKART[FFCellX - 1][FFCellY] > _kkart) && (IndexK <= Program.KKART[FFCellX - 1][FFCellY]) && (Program.CUTK[FFCellX - 1][FFCellY] > 0))
                {
                    velz += velz / MathF.Abs(velz) * 0.67F * Program.FloatMax(-Program.UK[FFCellX + 1][FFCellY][IndexK], 0);
                }
                //solid wall east of the particle
                if ((Program.KKART[FFCellX + 1][FFCellY] > _kkart) && (IndexK <= Program.KKART[FFCellX + 1][FFCellY]) && (Program.CUTK[FFCellX + 1][FFCellY] > 0))
                {
                    velz += velz / MathF.Abs(velz) * 0.67F * Program.FloatMax(Program.UK[FFCellX][FFCellY][IndexK], 0);
                }
                //solid wall south of the particle
                if ((Program.KKART[FFCellX][FFCellY - 1] > _kkart) && (IndexK <= Program.KKART[FFCellX][FFCellY - 1]) && (Program.CUTK[FFCellX][FFCellY - 1] > 0))
                {
                    velz += velz / MathF.Abs(velz) * 0.67F * Program.FloatMax(-Program.VK[FFCellX][FFCellY + 1][IndexK], 0);
                }
                //solid wall north of the particle
                if ((Program.KKART[FFCellX][FFCellY + 1] > _kkart) && (IndexK <= Program.KKART[FFCellX][FFCellY + 1]) && (Program.CUTK[FFCellX][FFCellY + 1] > 0))
                {
                    velz += velz / MathF.Abs(velz) * 0.67F * Program.FloatMax(Program.VK[FFCellX][FFCellY][IndexK], 0);
                }
                //******************************************************************************************************************** OETTL, 31 AUG 2016

                if (ObLength >= 0)
                {
                    Math.Clamp(velz, -3, 3);
                }
                else
                {
                    Math.Clamp(velz, -5, 5);
                }


                if (Program.GasBuoyancyMode != 0 && Program.DenseStates != null &&
                    nteil > 0 && nteil < Program.DenseStates.Length)
                {
                    float randSpread01 = (u_rg + 1) * RNG_Const;
                    float uxBefore = UXint;//20260308dbg
                    float uyBefore = UYint;//20260308dbg
                    float uzBefore = UZint;
                    ref Program.DenseGasState state = ref Program.DenseStates[nteil];
                    // update column bulk using current particle state at this idt
                    float rhoAirLocal = DenseRhoAir();
                    DenseColumnAccumulate(FFCellX, FFCellY, MathF.Max(0f, PartHeightAboveTerrain), state.rhoBulk, rhoAirLocal);
                    DenseColumnSample(FFCellX, FFCellY, rhoAirLocal, out float hLocal, out float rhoLocal, out float bLocal, out int cLocal);

                    state.HbulkLocal = hLocal;
                    state.RhoBulkLocal = rhoLocal;
                    state.BbulkLocal = bLocal;
                    state.ColumnCount = cLocal;

                    // mirror to dbg fields
                    state.dbgColHbulk = hLocal;
                    state.dbgColRhoBulk = rhoLocal;
                    state.dbgColBbulk = bLocal;
                    state.dbgColCount = cLocal;

                    ApplyDenseLightBuoyancy(ref UXint, ref UYint, ref UZint,
                                            ref state,
                                            idt, windge,
                                            PartHeightAboveTerrain,
                                            //FFCellX, FFCellY,
                                            randSpread01);
                    //20260305dbg
                    if (Program.GasBuoyancyMode == 1 &&
                        nteil <= 3 &&
                        timestep_number % 10 == 0)
                    {
                        float hLocalDbg = state.dbgColHbulk;
                        float hUsedDbg = MathF.Max(Program.DenseGroundEps, hLocalDbg > 0f ? hLocalDbg : state.cloudDepth);//20260308dbg

                        float UgsDbg = state.Ug * MathF.Cos(state.spreadAngle);//20260308dbg
                        float VgsDbg = state.Ug * MathF.Sin(state.spreadAngle);//20260308dbg

                        // 这一对能直接看出 dense 模块在水平上“实际加了多少”
                        float dUxDenseDbg = UXint - uxBefore;//20260308dbg
                        float dUyDenseDbg = UYint - uyBefore;//20260308dbg

                        // 调用后总水平速度（已包含 dense 修正）
                        float uTotDbg = MathF.Sqrt(Program.Pow2(UXint) + Program.Pow2(UYint));//20260308dbg
                        //DBG 追加 Udense/Ucarry//20260309
                        float uDenseDbg = MathF.Sqrt(UgsDbg * UgsDbg + VgsDbg * VgsDbg);
                        float uCarryDbg = MathF.Sqrt(state.up * state.up + state.vp * state.vp);

                        System.Console.WriteLine(
                            $"[DBG_DENSE] n={nteil} step={timestep_number} " +
                            $"h={PartHeightAboveTerrain:F3} H={state.dbgColHbulk:F3} " +
                            $"B={state.dbgColBbulk:E3} rhoB={state.dbgColRhoBulk:F4} cnt={state.dbgColCount} " +
                            $"phase={state.phase} hit={state.hitGroundOnce} " +
                            $"Ug={state.Ug:F3} UgPhys={state.dbgUgPhys:F3} UgCap={state.dbgUgCapByStep:F3} UgMax={Program.DenseUgMax:F3} " +
                            $"idt={state.dbgDt:F4} wind={state.dbgWind:F3} wp={state.dbgWpAbs:F3} " +
                            $"k={state.dbgK:E3} tau={state.dbgTau:F1}s relax={state.lastRelax:E3} " +
                            $"relaxChk={state.dbgRelaxCheck:E3} relax3s={state.dbgRelax3s:E3} wr={state.dbgWr:F3}" +
                            $"Hloc={hLocalDbg:F3} Huse={hUsedDbg:F3} " +
                            $"Ugs={UgsDbg:F3} Vgs={VgsDbg:F3} up={state.up:F3} vp={state.vp:F3} " +
                            $"dUx={dUxDenseDbg:F3} dUy={dUyDenseDbg:F3} Utot={uTotDbg:F3} " +
                            $"Udense={uDenseDbg:F3} Ucarry={uCarryDbg:F3} ");
                    }

                    if (timestep_number % 10000 == 0)
                    {
                        float uzAfter = UZint;
                        //Console.WriteLine($"[DBG] UZint delta={uzAfter - uzBefore:F6} before={uzBefore:F6} after={uzAfter:F6} velz={velz:F6} sum={velz + UZint:F6} dz={idt * (velz + UZint):F6}");
                    }
                }

                bool denseActive = Program.GasBuoyancyMode == 1 && Program.DenseStates != null &&
                   nteil > 0 && nteil < Program.DenseStates.Length &&
                   Program.DenseStates[nteil].phase != Program.DensePhase.Neutral;

                if (SourceType == Consts.SourceTypePoint && wpHurley > 0 && !useExplicitExitVelocityDirection) //20260521 1852 disable Hurley rise when explicit source direction is active
                {
                    //plume rise following Hurley, 2005 (TAPM)
                    //float epshurly = 1.5F * Program.Pow3(wpHurley) / (MeanHurley + 0.1F);
                    //float epsambiente = 0;
                    //epsambiente = eps;
                    DeltaZHurley = 0;
                    SigmaUpHurley = 0;

                    if (1.5F * Program.Pow3(wpHurley) / MeanHurley <= eps)
                    {
                        //if (nteil == 17110)//20251223 1416
                        //{
                            //Console.WriteLine($"[DBG] wpHurley_off n={nteil} step={timestep_number} reason=eps_exceeds plume_term wpHurley={wpHurley:F6} MeanHurley={MeanHurley:F6} eps={eps:F6}");//20251223 1323 鏃╂湡姝ユ暟鐭棩蹇楋紝瀹氫綅wpHurley=0
                        //}
                        wpHurley = 0;
                    }
                    else
                    {
                        float stab = 0;
                      
                        if (ObLength >= 0)
                        {
                            stab = 0.04F * MathF.Pow(2.73F, -ObLength * 0.05F);
                        }
                        
                        //plume-rise velocity
                        m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                        m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                        u_rg = (m_z << 16) + m_w;
                        u1_rg = (u_rg + 1) * RNG_Const;
                        m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                        m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                        u_rg = (m_z << 16) + m_w;
                        zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);
                        float windSpeedStandDev = 0;  
                        
                        if (auszeit < 3) // reduced wind fluctuation at the start of the plume rise for higher wind speeds
                        {
                            float maxFac = MathF.Max(8 - windge * 0.3F, 1.2F);
                            float minFac = MathF.Min(windge * 0.06F, 0.8F);
                            windSpeedStandDev = windge * Math.Clamp(1 + (windge * 0.35F + 0.25F) * zahl1, minFac, maxFac);
                        }
                        else
                        {
                            windSpeedStandDev = windge * Math.Clamp(1 + (windge * 0.31F + 0.25F) * zahl1, 0.1F, 8);
                        }
                                                
                        //Plume volume
                        GHurley += 2 * RHurley * (aHurley * Program.Pow2(wpHurley) + bHurley * windSpeedStandDev * wpHurley
                                                        + 0.1F * upHurley * MathF.Sqrt(0.5F * (Program.Pow2(velxold) + Program.Pow2(velyold)))) * idt;
                        FHurley += -0.035934F * stab * MHurley / upHurley * (0.4444444F * windSpeedStandDev + wpHurley) * idt;
                        //Momentum flux
                        MHurley += FHurley * idt;
                        //Plume radius
                        //RHurley = MathF.Sqrt((GHurley + FHurley / 9.8F) / upHurley);
                        RHurley = MathF.Sqrt((GHurley + FHurley * 0.1020408F) / upHurley);
                        upHurley = MathF.Sqrt(Program.Pow2(windSpeedStandDev) + Program.Pow2(wpHurley));
                        //standard deviations of velocity 
                        float sigmawpHurley = (aHurley * Program.Pow2(wpHurley) + bHurley * windSpeedStandDev * wpHurley) / (4.2F * upHurley);
                        SigmaUpHurley = 2 * sigmawpHurley;
                        float wpold = wpHurley;
                        wpHurley = Program.FloatMax(MHurley / GHurley, 0); // [m/s]
                        float wpmittel = (wpHurley + wpold) * 0.5F;

                        m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                        m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                        u_rg = (m_z << 16) + m_w;
                        zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);
                        
                        if (auszeit < 0.6F) // fix plume rise algorithm problems at the start of the plume rise
                        {
                            wpold *= MathF.Exp(-(auszeit + idt));
                            if (wpHurley < wpold) // [m/s] for "cold" stack sources at the start; wpold = exitVelocity  
                            {
                                if (MeanHurley < minVDI3782) // [m] plume rise height still below minimum plume rise height acc. to VDI 3782 equation 18
                                {
                                    wpHurley = MathF.Min(wpold, (minVDI3782 - MeanHurley) / idt); //[m/s] increase the plume rise height, until minimum plume rise height is reached
                                    wpmittel = wpHurley;
                                }
                            }
                            DeltaZHurley = MathF.Max(0, (wpmittel + sigmawpHurley * zahl1) * idt); // [m/s] * [s] = [m]
                        }
                        else
                        {
                            DeltaZHurley = Math.Max(0, (wpmittel + sigmawpHurley * zahl1) * idt); // [m/s] * [s] = [m]
                        }
                        MeanHurley += DeltaZHurley; // Height of plume delta_h(x) 
                    }

                    if (tunfak == Consts.ParticleIsNotAPortal)
                    {
                        float zBefore = zcoord_nteil;//20251223 0051
                        zcoord_nteil = zcoord_nteil + (velz + UZint) * idt + DeltaZHurley;
                    }
                    velzold = velz;
                }
                else
                {
                    if (tunfak == Consts.ParticleIsNotAPortal)
                    {
                        float zBefore = zcoord_nteil;//20251223 0051
                        zcoord_nteil += (velz + UZint) * idt;
                    }
                    velzold = velz;
                }

                if (Deposition_type > Consts.DepoOff && vsed > 0) // compute sedimentation for this particle
                {
                    float zBeforeSed = zcoord_nteil;
                    zcoord_nteil -= vsed * idt;
                }

                //control dispersion time of particles
                auszeit += idt;
                if (ISTATIONAER == Consts.TransientMode)
                {
                    if (auszeit >= tgesamt)
                    {
                        if (Program.IWET <= 2 && nteil == 1) // debug
                        {
                            Console.WriteLine($"[DBG] Conz5d write triggered: auszeit={auszeit:F3} tgesamt={tgesamt:F3}");
                        }

                        //particles from point sources with exit velocities are tracked until this velocity becomes almost zero
                        if (SourceType == Consts.SourceTypePoint)
                        {
                            if (wpHurley < 0.01)
                            {
                                Program.AccumulateIntraStepSnapshotFromConz5d(reflexion_flag, zcoord_nteil, AHint, masse, Area_cart, idt, xsi, eta, SG_nteil, auszeit);//20260402
                                TransientConcentration.Conz5dZeitschleife(reflexion_flag, zcoord_nteil, AHint, masse, Area_cart, idt, xsi, eta, SG_nteil);
                                //goto REMOVE_PARTICLE;
                                if (Program.ContinuousTraj) // keep particle for next outer step
                                {
                                    endOfStep = true;
                                    break; // end current step without removing particle
                                }
                                else
                                {
                                    goto REMOVE_PARTICLE;
                                }
                            }
                        }
                        //particles from portal sources with exit velocities are tracked until this velocity became zero
                        else if (SourceType == Consts.SourceTypePortal)
                        {
                            if (tunfak == Consts.ParticleIsNotAPortal)
                            {
                                Program.AccumulateIntraStepSnapshotFromConz5d(reflexion_flag, zcoord_nteil, AHint, masse, Area_cart, idt, xsi, eta, SG_nteil, auszeit);//20260402
                                TransientConcentration.Conz5dZeitschleife(reflexion_flag, zcoord_nteil, AHint, masse, Area_cart, idt, xsi, eta, SG_nteil);
                                //goto REMOVE_PARTICLE;
                                if (Program.ContinuousTraj) // keep particle for next outer step
                                {
                                    endOfStep = true;
                                    break; // end current step without removing particle
                                }
                                else
                                {
                                    goto REMOVE_PARTICLE;
                                }
                            }
                        }
                        else
                        {
                            Program.AccumulateIntraStepSnapshotFromConz5d(reflexion_flag, zcoord_nteil, AHint, masse, Area_cart, idt, xsi, eta, SG_nteil, auszeit);//20260402
                            TransientConcentration.Conz5dZeitschleife(reflexion_flag, zcoord_nteil, AHint, masse, Area_cart, idt, xsi, eta, SG_nteil);
                            //goto REMOVE_PARTICLE; // old behavior
                            if (Program.ContinuousTraj) // keep particle for next outer step
                            {
                                endOfStep = true;
                                break; // end current step without removing particle
                            }
                            else
                            {
                                goto REMOVE_PARTICLE;
                            }
                        }
                    }
                }

                //remove particles above maximum orography
                if (zcoord_nteil >= Program.ModelTopHeight)
                {
                    Program.Remove_Top++;//20260114
                    goto REMOVE_PARTICLE;
                }

                //new space coordinates due to turbulent movements
                float corx = velxold * idt;
                float cory = velyold * idt;

                /*
                 *    HORIZONTAL DIFFUSION ACCORDING TO ANFOSSI ET AL. (2006)
                 */

                //vertical interpolation of horizontal standard deviations of wind speed components
                (U0int, V0int) = IntStandCalculate(nteil, roughZ0, PartHeightAboveBuilding, windge, SigmaUpHurley);

                //determination of the meandering parameters according to Oettl et al. (2006)
                float param = 0;
                if (ISTATISTIK != Consts.MeteoSonic)
                {
                    param = Program.FloatMax(8.5F / Program.Pow2(windge + 1), 0F);
                }
                else
                {
                    param = Program.MWindMEander;
                }

                float T3 = 0;
                float termp = 0;
                float termq = 0;
                if (windge >= 2.5F || Program.MeanderingOff || param == 0 || Program.TAUS <= 600)
                {
                    T3 = 2 * Program.Pow2(V0int) / Program.C0x / eps;
                    termp = 1 / T3;
                }
                else
                {
                    termp = 1 / (2 * Program.Pow2(V0int) / Program.C0x / eps);
                    if (ISTATISTIK != Consts.MeteoSonic)
                    {
                        float Tstern = 200 * param + 350;
                        T3 = Tstern / 6.28F / (Program.Pow2(param) + 1) * param;
                    }
                    else if (ISTATISTIK == Consts.MeteoSonic)
                    {
                        T3 = Program.TWindMeander;
                    }

                    termq = param / (Program.Pow2(param) + 1) / T3;
                }

                //random numbers for horizontal turbulent velocities
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                u1_rg = (u_rg + 1) * RNG_Const;
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                float zuff1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

                //float zuff1 = (float)SimpleRNG.GetNormal();
                zuff1 = Math.Clamp(zuff1, -2f, 2f);

                //random numbers for horizontal turbulent velocities
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                u1_rg = (u_rg + 1) * RNG_Const;
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                u_rg = (m_z << 16) + m_w;
                float zuff2 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);
                zuff2 = Math.Clamp(zuff2, -2f, 2f);

                float velX = velxold + UXint;
                float velY = velyold + UYint;
                if (nteil <= 3 && Program.IWET == 1 && auszeit < 1.0f)//20260202
                {
                    System.Console.WriteLine(
                        $"[DBG] step n={nteil} x={xcoord_nteil:F3} y={ycoord_nteil:F3} z={zcoord_nteil:F3} " +
                        $"idt={idt:F4} UX={UXint:F3} UY={UYint:F3} UZ={UZint:F3} " +
                        $"vx={velxold:F3} vy={velyold:F3} vz={velz:F3}"
                    );
                }
                //horizontal Langevin equations
                {
                    float velx = velxold - (termp * velxold + termq * velyold) * idt + U0int * MathF.Sqrt(2 * termp * idt) * zuff1
                    + idt * (DUDX * velX + DUDY * velY + U0int * DSIGUDX)
                    + idt * (velxold / U0int * (DSIGUDX * velX + DSIGUDY * velY));

                    float vely = velyold - (termp * velyold - termq * velxold) * idt + V0int * MathF.Sqrt(2 * termp * idt) * zuff2
                    + idt * (DVDX * velX + DVDY * velY + V0int * DSIGVDY)
                    + idt * (velyold / V0int * (DSIGVDX * velX + DSIGVDY * velY));
                    velxold = velx;
                    velyold = vely;
                    if (float.IsNaN(velxold) || float.IsNaN(velyold) || float.IsInfinity(velxold) || float.IsInfinity(velyold))//20260114鏂板
                    {
                        System.Console.WriteLine($"[DBG] vel NaN/Inf n={nteil} U0int={U0int} V0int={V0int} termp={termp} termq={termq} idt={idt}");
                        goto REMOVE_PARTICLE;
                    }
                }
                //update coordinates
                if (tunfak == Consts.ParticleIsNotAPortal)
                {
                    float denseRandX = 0f;
                    float denseRandY = 0f;

                    if (Program.GasBuoyancyMode == 1 &&
                        Program.DenseStates != null &&
                        nteil > 0 &&
                        nteil < Program.DenseStates.Length)
                    {
                        ref Program.DenseGasState st = ref Program.DenseStates[nteil];
                        DenseHorzRandomStep(ref st, UXint, UYint, idt, ref m_z, ref m_w, out denseRandX, out denseRandY);
                    }

                    xcoord_nteil += idt * UXint + corx + denseRandX;
                    ycoord_nteil += idt * UYint + cory + denseRandY;
                    if (double.IsNaN(xcoord_nteil) || double.IsNaN(ycoord_nteil))//20260114鏂板
                    {
                        System.Console.WriteLine($"[DBG] xy NaN n={nteil} x={xcoord_nteil} y={ycoord_nteil} idt={idt} UXint={UXint} UYint={UYint} corx={corx} cory={cory}");
                        goto REMOVE_PARTICLE;
                    }
                }

                /*
                 *    TUNNEL MODULE
                 */
                if (tunfak == Consts.ParticleIsAPortal)
                {
                    //switch to standard dispersion
                    if ((Math.Abs(UXint - tunx) < 0.01) && (Math.Abs(UYint - tuny) < 0.01))
                    {
                        tunfak = Consts.ParticleIsNotAPortal;
                    }
                    if (tunpa <= 0)
                    {
                        tunfak = Consts.ParticleIsNotAPortal;
                    }

                    //particles within user-defined opposite lane are treated with standard dispersion (tunnel jet is destroyed there)
                    if (Program.TunnelOppLane == true)
                    {
                        if ((Program.OPP_LANE[FFCellX][FFCellY] == 1) && (PartHeightAboveTerrain <= 4))
                        {
                            tunfak = Consts.ParticleIsNotAPortal;

                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            u1_rg = (u_rg + 1) * RNG_Const;
                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

                            velxold = 2 * zahl1;

                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            u1_rg = (u_rg + 1) * RNG_Const;
                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

                            velyold = 2 * zahl1;

                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            u1_rg = (u_rg + 1) * RNG_Const;
                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

                            velzold = 2 * zahl1;

                            goto MOVE_FORWARD;
                        }
                    }

                    UXint = UXneu;
                    UYint = UYneu;

                    //maximum time step
                    idt = Program.GralDx / (tunpa + 0.05F) * 0.5F;
                    Math.Clamp(idt, 0.05F, 0.5F);

                    //idt = (float)Math.Min(Program.GralDx / (tunpa + 0.05F) * 0.5F, 0.5);
                    float tunpao = tunpa;

                    //dispersion time of each particle
                    auszeit += idt;

                    //calculation of cross- and along-wind
                    float windpa = (float)(UXint * tuncos + UYint * tunsin);
                    float windqu = (float)(-UXint * tunsin + UYint * tuncos);

                    //bending tunnel jet
                    if (windqu == 0)
                    {
                        windqu = 0.000001F;
                    }

                    float tunjettraeg = 1;
                    tunjettraeg = 5 * MathF.Exp(-0.005F * Program.TS_Area[Kenn_NTeil] * TS_ExitVelocity);
                    tunjettraeg = Program.FloatMax(tunjettraeg, 0.05F);
                    tunjettraeg = Math.Min(tunjettraeg, 1.5F);
                    tunqu = idt * 0.5F * Program.Pow2(windqu) * Math.Abs(windqu) / windqu * tunjettraeg;

                    //decelaration of the tunnel jet by friction/entrainment
                    tunjettraeg = 0.3F;
                    tunpa -= auszeit * idt * (tunpa - windpa) / (Program.Pow2(tunbreitalt)) * tunjettraeg;
                    if (tunpa < 0)
                    {
                        tunpa = 0;
                        tunfak = Consts.ParticleIsNotAPortal;
                    }

                    //calculation of jet stream velocities components
                    tunx = tunpa * tuncos - tunqu * tunsin;
                    tuny = tunpa * tunsin + tunqu * tuncos;

                    //coordinates of the central jet stream
                    xtun += idt * tunx;
                    ytun += idt * tuny;
                    tuncos = (float)(tunx / (Math.Sqrt(Program.Pow2(tunx) + Program.Pow2(tuny) + 0.001F)));
                    tunsin = (float)(tuny / (Math.Sqrt(Program.Pow2(tunx) + Program.Pow2(tuny) + 0.001F)));

                    //lateral widening of jet stream due to decelearation (mass continuity); but no narrowing when acceleration
                    if (tunpao > tunpa)
                    {
                        tunbreit = tunbreitalt * (Math.Abs(tunpao / (tunpa + 0.01F)));
                    }
                    yold = yneu;
                    yneu = yold + yold / tunbreitalt * (tunbreitalt - tunbreitalt);
                    tunbreitalt = Program.FloatMax(tunbreit, 0.000001F);

                    //linear increase of time-scale for jet-stream turbulence
                    float TLW = 2 * PartHeightAboveBuilding / (tunpa + 0.01F);

                    //time-scale for ambient-air turbulence
                    float sigmaz = 1.25F * Ustern;
                    float epsa = (Program.Pow3(Ustern) / 0.4F / Program.FloatMax(PartHeightAboveBuilding, Program.Z0Gramm[GrammCellX][GrammCellY]) *
                                         MathF.Pow(1 + 0.5F * MathF.Pow(PartHeightAboveBuilding / Math.Abs(ObLength), 0.8F), 1.8F));
                    float epst = 0.06F * Program.Pow2(tunpa) / TLW;

                    //in case of negative temperature differences (stable jet stream) the Lagrangian time-scale is reduced
                    if (TS_ExitTemperature < 0)
                    {
                        epst = Math.Max(epst / (10 * Math.Max(0.1F, Program.Pow2(TS_ExitTemperature))), epst * 0.01F);
                    }
                    else
                    {
                        epst *= (1 + TS_ExitTemperature * 0.5F) * 0.5F;
                    }

                    float TE = 2 * Program.Pow2(sigmaz) / Program.C0z / epsa;
                    float termw = 1 / TLW;

                    //temperature difference tunnel-jet <-> ambient air
                    m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                    m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                    u_rg = (m_z << 16) + m_w;
                    u1_rg = (u_rg + 1) * RNG_Const;
                    m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                    m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                    u_rg = (m_z << 16) + m_w;
                    zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

                    float zuff4 = idt * zahl1;
                    buoy += -buoy * termw * idt + zuff4 * MathF.Sqrt(epst);

                    //new z-coordinate of particles
                    zturb = buoy * idt;
                    if (TS_ExitTemperature >= 0)
                    {
                        zturb = buoy * idt + velz * idt;
                    }

                    zcoord_nteil += zturb;

                    //mass-continuity in bending areas
                    double xoldy = xneuy;
                    double yoldy = yneuy;
                    xneuy = xtun - yteilchen * tunsin;
                    yneuy = ytun + yteilchen * tuncos;

                    //particle coordinates

                    xturb += corx;
                    yturb += cory;
                    xcoord_nteil = xtun - yneu * tunsin + xturb;
                    ycoord_nteil = ytun + yneu * tuncos + yturb;

                    //coordinate transformation towards model grid
                    xsi = xcoord_nteil - IKOOAGRAL;
                    eta = ycoord_nteil - JKOOAGRAL;
                    xsi1 = xcoord_nteil - Program.GrammWest;
                    eta1 = ycoord_nteil - Program.GrammSouth;
                    if ((eta <= Program.EtaMinGral) || (xsi <= Program.XsiMinGral) || (eta >= Program.EtaMaxGral) || (xsi >= Program.XsiMaxGral))
                    {
                        if (Program.Remove_OutDomain < 10)
                        {
                            string outreason =
                                (xsi <= Program.XsiMinGral) ? "xsi<=min" :
                                (xsi >= Program.XsiMaxGral) ? "xsi>=max" :
                                (eta <= Program.EtaMinGral) ? "eta<=min" :
                                "eta>=max";

                            System.Console.WriteLine(
                                $"[DBG] OUT[{outreason}] n={nteil} IWET={Program.IWET} " +
                                $"x={xcoord_nteil} y={ycoord_nteil} z={zcoord_nteil} " +
                                $"PartH={(zcoord_nteil - AHint)} AHint={AHint} " +
                                $"xsi={xsi} eta={eta} XsiMin={Program.XsiMinGral} XsiMax={Program.XsiMaxGral} " +
                                $"EtaMin={Program.EtaMinGral} EtaMax={Program.EtaMaxGral}"
                            );
                        }  
                        Program.Remove_OutDomain++;//20260114
                        goto REMOVE_PARTICLE;
                    }

                    int IndexIalt = 1;
                    int IndexJalt = 1;
                    if (topo == Consts.TerrainAvailable)
                    {
                        //with topography
                        IndexIalt = FFCellX;
                        IndexJalt = FFCellY;
                        //index in the flow field grid
                        FFCellX = (int)(xsi * FFGridXRez) + 1;
                        FFCellY = (int)(eta * FFGridYRez) + 1;
                        if ((FFCellX < 1) || (FFCellX > Program.NII) || (FFCellY < 1) || (FFCellY > Program.NJJ))
                        {
                            goto REMOVE_PARTICLE;
                        }

                        //index in the GRAMM grid
                        GrammCellX = Math.Clamp((int)(xsi1 * GrammGridXRez) + 1, 1, Program.NX);
                        GrammCellY = Math.Clamp((int)(eta1 * GrammGridYRez) + 1, 1, Program.NY);

                        //interpolated orography
                        AHintold = AHint;
                        AHint = Program.AHK[FFCellX][FFCellY];

                        PartHeightAboveTerrain = zcoord_nteil - AHint;
                        PartHeightAboveBuilding = PartHeightAboveTerrain;
                    }
                    else
                    {
                        //flat terrain
                        FFCellX = (int)(xsi * FFGridXRez) + 1;
                        FFCellY = (int)(eta * FFGridYRez) + 1;
                        //if ((FFCellX > Program.NII) || (FFCellX > Program.NJJ) || (FFCellY < 1) || (FFCellY < 1))//鍘熷鐗堟湰
                        if ((FFCellX > Program.NII) || (FFCellY > Program.NJJ) || (FFCellX < 1) || (FFCellY < 1))//20260114
                        {
                            if (Program.Remove_FFCell < 5)//20260114
                            {
                                System.Console.WriteLine($"[DBG] FFCell drop x={xcoord_nteil:F1} y={ycoord_nteil:F1} FFCellX={FFCellX} FFCellY={FFCellY} NII={Program.NII} NJJ={Program.NJJ}");
                            }
                            Program.Remove_FFCell++; // 涓存椂璁℃暟-20260114
                            goto REMOVE_PARTICLE;
                        }

                        //interpolated orography
                        AHintold = AHint;
                        PartHeightAboveTerrain = zcoord_nteil;
                        PartHeightAboveBuilding = PartHeightAboveTerrain;
                    }

                    //remove particles above maximum orography
                    if (zcoord_nteil >= Program.ModelTopHeight)
                    {
                        Program.Remove_Top++;//20260114
                        goto REMOVE_PARTICLE;
                    }

                    //in case that orography is taken into account and if there is no particle, then, the particle trajectories are following the terrain
                    if (topo == Consts.TerrainAvailable)
                    {
                        // no terrain following trajectory, if the recent or the previous cell is a building!
                        if ((Program.CUTK[FFCellX][FFCellY] + Program.CUTK[FFCellXPrev][FFCellYPrev]) == 0 )
                        {
                            zcoord_nteil += (AHint - AHintold);
                        }
                    }

                    //time-step correction to ensure mass-conservation for bending tunnel jets
                    geschwneu = (float)(Math.Sqrt(Program.Pow2(xneuy - xoldy) + Program.Pow2(yneuy - yoldy)) / idt);
                    if (geschwalt == 0)
                    {
                        geschwalt = geschwneu;
                    }

                    faktzeit = Program.FloatMax(faktzeit * geschwneu / (geschwalt + 0.001F), 0.001F);
                    idt *= faktzeit;
                    if (idt < 0.01)
                    {
                        idt = 0.01f; // 29.12.2017 Kuntner: avoid extremely low time-steps at extreme wind speeds
                    }

                    geschwalt = geschwneu;

                    //remove particles above GRAMM orography
                    if (topo == Consts.TerrainAvailable)
                    {
                        if (zcoord_nteil >= Program.ModelTopHeight)
                        {
                            Program.Remove_Top++;//20260114
                            goto REMOVE_PARTICLE;
                        }
                        //if (zcoord_nteil >= Program.AHMAX + 290)
                        //	goto REMOVE_PARTICLE;
                    }

                    //reflexion of particles from tunnel portals
                    reflexion_flag = Consts.ParticleNotReflected;
                    int IndexKOld = 0;
                    if (topo == Consts.TerrainAvailable)
                    {
                        if (Program.BuildingsExist == true)
                        {
                            IndexKOld = IndexK;
                            IndexK = BinarySearch(zcoord_nteil - Program.AHMIN); //19.05.25 Ku

                            if ((IndexK <= Program.KKART[FFCellX][FFCellY]) && (Program.CUTK[FFCellX][FFCellY] > 0))
                            {
                                // compute deposition according to VDI 3945 for this particle- add deposition to Depo_conz[][][]
                                if (Deposition_type > Consts.DepoOff && depo_reflection_counter >= 0)
                                {
                                    float Pd1 = Pd(varw, vsed, vdep, Deposition_type, FFCellX, FFCellY);
                                    int ik = (int)(xsi * ConcGridXRez) + 1;
                                    int jk = (int)(eta * ConcGridYRez) + 1;
                                    double[] depo_L = Program.Depo_conz[ik][jk];
                                    double conc = masse * Pd1 * area_rez_fac;
                                    lock (depo_L)
                                    {
                                        depo_L[SG_nteil] += conc;
                                        depo_events++;//20251223 娌夐檷绱姞
                                    }
                                    masse -= masse * Pd1;
                                    depo_reflection_counter = -2; // block deposition 1 for 1 timestep
                                    
                                    if (Program.ContinuousTraj)
                                    {
                                        if (masse < 1e-12) masse = 1e-12; // continuous mode: keep tiny residual mass
                                    }
                                    else
                                    {
                                        if (masse <= 0)
                                        {
                                            Program.Remove_MassZero++;//20260114
                                            goto REMOVE_PARTICLE;
                                        }
                                    }
                                }

                                tunfak = Consts.ParticleIsNotAPortal;
                                reflexion_flag = Consts.ParticleReflected;

                                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                                u_rg = (m_z << 16) + m_w;
                                u1_rg = (u_rg + 1) * RNG_Const;
                                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                                u_rg = (m_z << 16) + m_w;
                                zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

                                velzold = (MathF.Sqrt(varw) + 0.82F * tunpa) * zahl1;

                                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                                u_rg = (m_z << 16) + m_w;
                                u1_rg = (u_rg + 1) * RNG_Const;
                                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                                u_rg = (m_z << 16) + m_w;
                                zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);
                                velxold = (U0int + 0.82F * tunpa) * zahl1;

                                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                                u_rg = (m_z << 16) + m_w;
                                u1_rg = (u_rg + 1) * RNG_Const;
                                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                                u_rg = (m_z << 16) + m_w;
                                zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);
                                velyold = (V0int + 0.82F * tunpa) * zahl1;

                                tunpa = 0.01F;
                                goto MOVE_FORWARD;
                            }
                        }
                    }
                    else if ((topo == Consts.TerrainFlat) && (Program.BuildingsExist == true))
                    {
                        //flat terrain with buildings
                        IndexKOld = IndexK;
                        IndexK = BinarySearch(zcoord_nteil); //19.05.25 Ku

                        if (IndexK <= Program.KKART[FFCellX][FFCellY])
                        {
                            // compute deposition according to VDI 3945 for this particle- add deposition to Depo_conz[][][]
                            if (Deposition_type > Consts.DepoOff && depo_reflection_counter >= 0)
                            {
                                float Pd1 = Pd(varw, vsed, vdep, Deposition_type, FFCellX, FFCellY);
                                int ik = (int)(xsi * ConcGridXRez) + 1;
                                int jk = (int)(eta * ConcGridYRez) + 1;
                                double[] depo_L = Program.Depo_conz[ik][jk];
                                double conc = masse * Pd1 * area_rez_fac;
                                lock (depo_L)
                                {
                                    depo_L[SG_nteil] += conc;
                                    depo_events++;//20251223 绱姞
                                }
                                masse -= masse * Pd1;
                                depo_reflection_counter = -2; // block deposition 1 for 1 timestep

                                if (Program.ContinuousTraj)
                                {
                                    if (masse < 1e-12) masse = 1e-12; // continuous mode: keep tiny residual mass
                                }
                                else
                                {
                                    if (masse <= 0)
                                    {
                                        Program.Remove_MassZero++;//20260114
                                        goto REMOVE_PARTICLE;
                                    }
                                }
                            }

                            tunfak = Consts.ParticleIsNotAPortal;
                            reflexion_flag = Consts.ParticleReflected;

                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            u1_rg = (u_rg + 1) * RNG_Const;
                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);
                            velzold = (MathF.Sqrt(varw) + 0.82F * tunpa) * zahl1;


                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            u1_rg = (u_rg + 1) * RNG_Const;
                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

                            velxold = (U0int + 0.82F * tunpa) * zahl1;

                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            u1_rg = (u_rg + 1) * RNG_Const;
                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

                            velyold = (V0int + 0.82F * tunpa) * zahl1;

                            tunpa = 0.01F;
                            goto MOVE_FORWARD;
                        }
                    }

                    //reflexion at the surface and at the top of the boundary layer
                    if (PartHeightAboveTerrain <= 0)
                    {
                        zcoord_nteil = AHint - PartHeightAboveTerrain + 0.01F;
                        PartHeightAboveTerrain = zcoord_nteil - AHint;
                        PartHeightAboveBuilding = PartHeightAboveTerrain;
                        
                        // compute deposition according to VDI 3945 for this particle- add deposition to Depo_conz[][][]
                        if (Deposition_type > Consts.DepoOff && depo_reflection_counter >= 0)
                        {
                            float Pd1 = Pd(varw, vsed, vdep, Deposition_type, FFCellX, FFCellY);
                            int ik = (int)(xsi * ConcGridXRez) + 1;
                            int jk = (int)(eta * ConcGridYRez) + 1;
                            double[] depo_L = Program.Depo_conz[ik][jk];
                            double conc = masse * Pd1 * area_rez_fac;
                            lock (depo_L)
                            {
                                depo_L[SG_nteil] += conc;
                                depo_events++;//20251223
                            }
                            masse -= masse * Pd1;
                            depo_reflection_counter = -2; // block deposition 1 for 1 timestep

                            if (Program.ContinuousTraj)
                            {
                                if (masse < 1e-12) masse = 1e-12; // continuous mode: keep tiny residual mass
                            }
                            else
                            {
                                if (masse <= 0)
                                {
                                    Program.Remove_MassZero++;//20260114
                                    goto REMOVE_PARTICLE;
                                }
                            }
                        }
                    }

                    if (PartHeightAboveBuilding >= blh && Program.ISTATIONAER != Consts.TransientMode) //26042020 (Ku): removed for transient mode -> particles shoulb be tracked above blh 
                    {
                        goto REMOVE_PARTICLE;
                    }

                    goto MOVE_TO_CONCENTRATIONCALCULATION;

                } //END OF TUNNEL MODUL
                
                //in case that particle was reflected, flag is set to 1 and subsequently concentrations are not computed
                reflexion_flag = Consts.ParticleNotReflected;
                int back = 1;
                while (back == 1)
                {
                    //coordinate transformation towards model grid
                    xsi = xcoord_nteil - IKOOAGRAL;
                    eta = ycoord_nteil - JKOOAGRAL;
                    xsi1 = xcoord_nteil - Program.GrammWest;
                    eta1 = ycoord_nteil - Program.GrammSouth;
                    if ((eta <= Program.EtaMinGral) || (xsi <= Program.XsiMinGral) || (eta >= Program.EtaMaxGral) || (xsi >= Program.XsiMaxGral))
                    {
                        if (Program.Remove_OutDomain < 10)
                        {
                            string outreason =
                                (xsi <= Program.XsiMinGral) ? "xsi<=min" :
                                (xsi >= Program.XsiMaxGral) ? "xsi>=max" :
                                (eta <= Program.EtaMinGral) ? "eta<=min" :
                                "eta>=max";

                            System.Console.WriteLine(
                                $"[DBG] OUT[{outreason}] n={nteil} IWET={Program.IWET} " +
                                $"x={xcoord_nteil} y={ycoord_nteil} z={zcoord_nteil} " +
                                $"PartH={(zcoord_nteil - AHint)} AHint={AHint} " +
                                $"xsi={xsi} eta={eta} XsiMin={Program.XsiMinGral} XsiMax={Program.XsiMaxGral} " +
                                $"EtaMin={Program.EtaMinGral} EtaMax={Program.EtaMaxGral}"
                            );
                        }     
                        Program.Remove_OutDomain++;//20260114
                        goto REMOVE_PARTICLE;
                    }

                    if (topo == Consts.TerrainAvailable)
                    {
                        //with topography
                        FFCellXPrev = FFCellX;
                        FFCellYPrev = FFCellY;
                        FFCellX = (int)(xsi * FFGridXRez) + 1;
                        FFCellY = (int)(eta * FFGridYRez) + 1;
                        GrammCellX = Math.Clamp((int)(xsi1 * GrammGridXRez) + 1, 1, Program.NX);
                        GrammCellY = Math.Clamp((int)(eta1 * GrammGridYRez) + 1, 1, Program.NY);

                        if ((FFCellX > Program.NII) || (FFCellY > Program.NJJ) || (FFCellX < 1) || (FFCellY < 1))
                        {
                            goto REMOVE_PARTICLE;
                        }
                        //19.05.25 Ku: allow a relative step along the terrain near the start point one times
                        if (auszeit < 20 && TerrainStepAllowed && SourceType > 1) // Near the source, terrain correction 1 times, line and area sources only
                        {
                            int IndexXDelta = 2 * FFCellX - FFCellXPrev;
                            int IndexYDelta = 2 * FFCellY - FFCellYPrev;
                            if (IndexXDelta > 1 && IndexXDelta < Program.NII && IndexYDelta > 1 && IndexYDelta < Program.NJJ)
                            {
                                if (Program.CUTK[FFCellX][FFCellY] == 0 && Program.CUTK[FFCellXPrev][FFCellYPrev] == 0 && Program.CUTK[IndexXDelta][IndexYDelta] == 0) //No buildings
                                {
                                    int KKARTOld = Program.KKART[FFCellXPrev][FFCellYPrev];
                                    int KKARTNew = Program.KKART[IndexXDelta][IndexYDelta];

                                    if (KKARTOld != KKARTNew &&                 //Terrain step
                                        Math.Abs(KKARTOld - KKARTNew) < 3)     // step < 3
                                    {
                                        float AHK = Program.AHK[IndexXDelta][IndexYDelta];
                                        float AHKold = Program.AHK[FFCellXPrev][FFCellYPrev];

                                        if (Math.Abs(AHK - AHKold) > Program.GralDx * 0.5F) // vertical step > Concentration raster / 2
                                        {
                                            float faktor = Math.Abs(zcoord_nteil - Program.AHK[FFCellX][FFCellY]) / Program.DZK[KKARTNew];
                                            if (faktor < 3) // avoid overflow and increase perf.
                                            {
                                                int KKARTRecent = Program.KKART[FFCellX][FFCellY];
                                                faktor = MathF.Exp(-faktor * faktor * faktor * faktor * 0.01F);
                                                if (KKARTRecent == KKARTOld && KKARTNew < KKARTRecent) // particle is still at a higher cell
                                                {
                                                    faktor = 0F;
                                                }
                                                else if (KKARTNew < KKARTOld) // reduce correction at a downstep
                                                {
                                                    faktor *= 0.5F;
                                                }

                                                if (faktor > 0)
                                                {
                                                    zcoord_nteil += (AHK - Program.AHK[FFCellXPrev][FFCellYPrev]) * faktor * 0.8F;
                                                    if (zcoord_nteil < Program.AHK[FFCellX][FFCellY])
                                                    {
                                                        zcoord_nteil = Program.AHK[FFCellX][FFCellY] + 0.2F;
                                                    }
                                                    TerrainStepAllowed = false;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        //flat terrain
                        FFCellXPrev = FFCellX;
                        FFCellYPrev = FFCellY;
                        FFCellX = (int)(xsi * FFGridXRez) + 1;
                        FFCellY = (int)(eta * FFGridYRez) + 1;
                        if ((FFCellX > Program.NII) || (FFCellY > Program.NJJ) || (FFCellX < 1) || (FFCellY < 1))
                        {
                            if (Program.Remove_FFCell < 5)//20260114
                            {
                                System.Console.WriteLine($"[DBG] FFCell drop x={xcoord_nteil:F1} y={ycoord_nteil:F1} FFCellX={FFCellX} FFCellY={FFCellY} NII={Program.NII} NJJ={Program.NJJ}");
                            }
                            Program.Remove_FFCell++; // 涓存椂璁℃暟-20260114
                            goto REMOVE_PARTICLE;
                        }
                    }

                    //reflexion of particles within buildings and at the surface
                    back = 0;
                    if (Program.BuildingsExist == true || topo == Consts.TerrainAvailable)
                    {
                        int IndexKOld = IndexK;
                        IndexK = BinarySearch(zcoord_nteil - Program.AHMIN); //19.05.25 Ku
                        
                        // Particle below building or terrain
                        if (IndexK <= Program.KKART[FFCellX][FFCellY])
                        {
                            // The particle enters a building or terrain step horizontally
                            if ((IndexK > Program.KKART[FFCellXPrev][FFCellYPrev]) && (IndexK == IndexKOld))
                            {
                                if (idt * UXint + corx <= 0)
                                {
                                    xcoord_nteil = 2 * (IKOOAGRAL + FFCellX * FFGridX) - xcoord_nteil + 0.01F;
                                }
                                else
                                {
                                    xcoord_nteil = 2 * (IKOOAGRAL + (FFCellX - 1) * FFGridX) - xcoord_nteil - 0.01F;
                                }

                                if (idt * UYint + cory <= 0)
                                {
                                    ycoord_nteil = 2 * (JKOOAGRAL + FFCellY * FFGridY) - ycoord_nteil + 0.01F;
                                }
                                else
                                {
                                    ycoord_nteil = 2 * (JKOOAGRAL + (FFCellY - 1) * FFGridY) - ycoord_nteil - 0.01F;
                                }
                            }
                            // Reflect the particle in x direction
                            else if ((IndexK > Program.KKART[FFCellXPrev][FFCellY]) && (IndexK == IndexKOld))
                            {
                                if (idt * UXint + corx <= 0)
                                {
                                    xcoord_nteil = 2 * (IKOOAGRAL + FFCellX * FFGridX) - xcoord_nteil + 0.01F;
                                }
                                else
                                {
                                    xcoord_nteil = 2 * (IKOOAGRAL + (FFCellX - 1) * FFGridX) - xcoord_nteil - 0.01F;
                                }
                            }
                            // Reflect the particle in y direction
                            else if ((IndexK > Program.KKART[FFCellX][FFCellYPrev]) && (IndexK == IndexKOld))
                            {
                                if (idt * UYint + cory <= 0)
                                {
                                    ycoord_nteil = 2 * (JKOOAGRAL + FFCellY * FFGridY) - ycoord_nteil + 0.01F;
                                }
                                else
                                {
                                    ycoord_nteil = 2 * (JKOOAGRAL + (FFCellY - 1) * FFGridY) - ycoord_nteil - 0.01F;
                                }
                            }
                            // Particle comes from z direction on the same cell 
                            else if ((IndexKOld > Program.KKART[FFCellX][FFCellY]) && (FFCellX == FFCellXPrev) && (FFCellY == FFCellYPrev))
                            {
                                zcoord_nteil = 2 * Program.AHK[FFCellX][FFCellY] - zcoord_nteil + 0.01F;
                                PartHeightAboveTerrain = zcoord_nteil - AHint;
                                if (topo == Consts.TerrainAvailable)
                                {
                                    PartHeightAboveBuilding = PartHeightAboveTerrain;
                                }
                                // velzold = -velzold; // replaced by dense reflection
                                if (Program.GasBuoyancyMode == 1 && Program.DenseStates != null &&
                                    nteil > 0 && nteil < Program.DenseStates.Length)
                                {
                                    ref Program.DenseGasState st = ref Program.DenseStates[nteil];
                                    ApplyDenseBottomReflection(ref velzold, ref st);
                                }
                                else
                                {
                                    velzold = -velzold;
                                }

                                
                                // compute deposition according to VDI 3945 for this particle- add deposition to Depo_conz[][][]
                                if (Deposition_type > Consts.DepoOff && depo_reflection_counter >= 0)
                                {
                                    float Pd1 = Pd(varw, vsed, vdep, Deposition_type, FFCellX, FFCellY);
                                    int ik = (int)(xsi * ConcGridXRez) + 1;
                                    int jk = (int)(eta * ConcGridYRez) + 1;
                                    double[] depo_L = Program.Depo_conz[ik][jk];
                                    double conc = masse * Pd1 * area_rez_fac;
                                    lock (depo_L)
                                    {
                                        depo_L[SG_nteil] += conc;
                                        depo_events++;//20251223
                                    }
                                    masse -= masse * Pd1;
                                    depo_reflection_counter = -2; // block deposition 1 for the entire reflection algorithm and 1 timestep

                                    if (Program.ContinuousTraj)
                                    {
                                        if (masse < 1e-12) masse = 1e-12; // continuous mode: keep tiny residual mass
                                    }
                                    else
                                    {
                                        if (masse <= 0)
                                        {
                                            Program.Remove_MassZero++;//20260114
                                            goto REMOVE_PARTICLE;
                                        }
                                    }
                                }
                            }
                            // Particle comes from the z direction (above the recent cell) and from another cell
                            else if ((IndexKOld > Program.KKART[FFCellX][FFCellY]))
                            {
                                if (Deposition_type > Consts.DepoOff && depo_reflection_counter >= 0)
                                {
                                    float Pd1 = Pd(varw, vsed, vdep, Deposition_type, FFCellX, FFCellY);
                                    int ik = (int)(xsi * ConcGridXRez) + 1;
                                    int jk = (int)(eta * ConcGridYRez) + 1;
                                    double[] depo_L = Program.Depo_conz[ik][jk];
                                    double conc = masse * Pd1 * area_rez_fac;
                                    lock (depo_L)
                                    {
                                        depo_L[SG_nteil] += conc;
                                        depo_events++;//20251223
                                    }
                                    masse -= masse * Pd1;
                                    depo_reflection_counter = -2; // block deposition 1 for the entire reflection algorithm and 1 timestep 

                                    if (Program.ContinuousTraj)
                                    {
                                        if (masse < 1e-12) masse = 1e-12; // continuous mode: keep tiny residual mass
                                    }
                                    else
                                    {
                                        if (masse <= 0)
                                        {
                                            Program.Remove_MassZero++;//20260114
                                            goto REMOVE_PARTICLE;
                                        }
                                    }
                                }
                                float terrainHeight = Program.AHK[FFCellX][FFCellY];
                                zcoord_nteil = MathF.Max(terrainHeight, 2 * terrainHeight - zcoord_nteil + 0.01F);
                            }
                            // Particle comes from the z direction (below the recent cell, particle direction up or down) and from another cell
                            else
                            {
                                // compute deposition according to VDI 3945 for this particle- add deposition to Depo_conz[][][]
                                if (Deposition_type > Consts.DepoOff && depo_reflection_counter >= 0)
                                {
                                    float Pd1 = Pd(varw, vsed, vdep, Deposition_type, FFCellX, FFCellY);
                                    int ik = (int)(xsi * ConcGridXRez) + 1;
                                    int jk = (int)(eta * ConcGridYRez) + 1;
                                    double[] depo_L = Program.Depo_conz[ik][jk];
                                    double conc = masse * Pd1 * area_rez_fac;
                                    lock (depo_L)
                                    {
                                        depo_L[SG_nteil] += conc;
                                        depo_events++;//20251223
                                    }
                                    masse -= masse * Pd1;
                                    depo_reflection_counter = -2; // block deposition 1 for the entire reflection algorithm and 1 timestep

                                    if (Program.ContinuousTraj)
                                    {
                                        if (masse < 1e-12) masse = 1e-12; // continuous mode: keep tiny residual mass
                                    }
                                    else
                                    {
                                        if (masse <= 0)
                                        {
                                            Program.Remove_MassZero++;//20260114
                                            goto REMOVE_PARTICLE;
                                        }
                                    }
                                }

                                if (tunpa == 0)
                                {
                                    // Terrain following particles at small steps if there is no building
                                    if (Program.CUTK[FFCellX][FFCellY] < 1 && (Program.AHK[FFCellX][FFCellY] - zcoord_nteil) < 10)
                                    {
                                        // Move particle above the new terrain surface
                                        zcoord_nteil = Program.AHK[FFCellX][FFCellY] + MathF.Abs((idt * (UZint + velz) + DeltaZHurley));
                                    }
                                    else
                                    // Below Building or high terrain step
                                    {
                                        // reset horizontal coordinates
                                        double deltax = xcoord_nteil - xcoord_nteil_Prev;
                                        double deltay = ycoord_nteil - ycoord_nteil_Prev;
                                        
                                        xcoord_nteil -= deltax * 1.05;
                                        ycoord_nteil -= deltay * 1.05;
                                        //zcoord_nteil -= (idt * (UZint + velz) + DeltaZHurley);

                                        FFCellXPrev = FFCellX;
                                        FFCellYPrev = FFCellY;

                                        FFCellX = (int)((xcoord_nteil - IKOOAGRAL) * FFGridXRez) + 1;
                                        FFCellY = (int)((ycoord_nteil - JKOOAGRAL) * FFGridYRez) + 1;
                                        if ((FFCellX > Program.NII) || (FFCellY > Program.NJJ) || (FFCellX < 1) || (FFCellY < 1))
                                        {
                                            goto REMOVE_PARTICLE;
                                        }
                                        float newTerrainHeight = Program.AHK[FFCellX][FFCellY];
                                        // if below new terrain -> multiple reflections
                                        if (zcoord_nteil < newTerrainHeight)
                                        {
                                            zcoord_nteil = 2 * newTerrainHeight - zcoord_nteil;
                                        }
                                        if (topo == Consts.TerrainAvailable)
                                        {
                                            PartHeightAboveTerrain = zcoord_nteil - newTerrainHeight;
                                            PartHeightAboveBuilding = PartHeightAboveTerrain;
                                        }
                                        else
                                        {
                                            PartHeightAboveTerrain = zcoord_nteil;
                                            PartHeightAboveBuilding = zcoord_nteil;
                                        }                                  
                                    }
                                }

                                // Change direction of the vertical velocity
                                int vorzeichen = -1;
                                if (velzold < 0)
                                {
                                    vorzeichen = 1;
                                }

                                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                                u_rg = (m_z << 16) + m_w;
                                u1_rg = (u_rg + 1) * RNG_Const;
                                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                                u_rg = (m_z << 16) + m_w;
                                zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

                                velzold = MathF.Abs(zahl1 * MathF.Sqrt(varw)) * vorzeichen;
                            }

                            // Change direction of the horizontal particle velocity
                            int vorzeichen1 = -1;
                            if (velxold < 0)
                            {
                                vorzeichen1 = 1;
                            }

                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            u1_rg = (u_rg + 1) * RNG_Const;
                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);
                            
                            velxold = MathF.Abs(zahl1 * U0int * 3) * vorzeichen1;
                            
                            vorzeichen1 = -1;
                            if (velyold < 0)
                            {
                                vorzeichen1 = 1;
                            }

                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            u1_rg = (u_rg + 1) * RNG_Const;
                            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                            u_rg = (m_z << 16) + m_w;
                            zahl1 = MathF.Sqrt(-2F * MathF.Log(u1_rg)) * MathF.Sin(Pi2F * (u_rg + 1) * RNG_Const);

                            velyold = MathF.Abs(zahl1 * V0int * 3) * vorzeichen1;
                            
                            back = 1;
                            idt = Program.FloatMax(idt * 0.5F, 0.05F);
                            if (tunpa > 0)
                            {
                                goto MOVE_FORWARD;
                            }
                        }
                    }

                    if (back > 0)
                    {
                        reflexion_flag = Consts.ParticleReflected;
                        reflexion_number++;
                        //in rare cases, particles get trapped near solid boundaries. Such particles are abandoned after a certain number of reflexions
                        if (reflexion_number > Max_Reflections)
                        {
                            goto REMOVE_PARTICLE;
                        }
                    }
                }   //END OF THE REFLEXION ALGORITHM
                depo_reflection_counter = 0;

                //interpolation of orography
                if (topo == Consts.TerrainAvailable)
                {
                    if ((FFCellX < 1) || (FFCellX > Program.NII) || (FFCellY < 1) || (FFCellY > Program.NJJ))
                    {
                        goto REMOVE_PARTICLE;
                    }

                    AHint = Program.AHK[FFCellX][FFCellY];

                    PartHeightAboveTerrain = zcoord_nteil - AHint;
                    PartHeightAboveBuilding = PartHeightAboveTerrain;
                }
                else
                {
                    if ((FFCellX < 1) || (FFCellX > Program.NII) || (FFCellY < 1) || (FFCellY > Program.NJJ))
                    {
                        goto REMOVE_PARTICLE;
                    }
                    AHint = 0;
                    PartHeightAboveTerrain = zcoord_nteil;
                    PartHeightAboveBuilding = PartHeightAboveTerrain;
                }

                //reflexion at the surface and at the top of the boundary layer               
                if (Program.AdaptiveRoughnessMax > 0)
                {
                    ObL = Program.OLGral[FFCellX][FFCellY];
                }
                else
                {
                    ObL = Program.Ob[GrammCellX][GrammCellY];
                }

                if ((PartHeightAboveBuilding >= blh) && (ObL < 0))
                {
                    if ((topo == Consts.TerrainAvailable) && (UZint >= 0) && Program.ISTATIONAER != Consts.TransientMode) //26042020 (Ku): removed for transient mode -> particles shoulb be tracked above blh
                    {
                        goto REMOVE_PARTICLE;
                    }

                    zcoord_nteil = blh + AHint - (PartHeightAboveBuilding - blh) - 0.01F;
                    velzold = -velzold;
                    PartHeightAboveTerrain = zcoord_nteil - AHint;
                    PartHeightAboveBuilding = PartHeightAboveTerrain;
                }

                if (PartHeightAboveTerrain <= 0)
                {
                    zcoord_nteil = AHint - PartHeightAboveTerrain + 0.01F;
                    // velzold = -velzold; // replaced by dense reflection
                    if (Program.GasBuoyancyMode == 1 && Program.DenseStates != null &&
                        nteil > 0 && nteil < Program.DenseStates.Length)
                    {
                        ref Program.DenseGasState st = ref Program.DenseStates[nteil];
                        ApplyDenseBottomReflection(ref velzold, ref st);
                    }
                    else
                    {
                        velzold = -velzold;
                    }
                    PartHeightAboveTerrain = zcoord_nteil - AHint;
                    PartHeightAboveBuilding = PartHeightAboveTerrain;
                   
                    // compute deposition according to VDI 3945 for this particle- add deposition to Depo_conz[][][]
                    if (Deposition_type > Consts.DepoOff && depo_reflection_counter >= 0)
                    {
                        float Pd1 = Pd(varw, vsed, vdep, Deposition_type, FFCellX, FFCellY);
                        int ik = (int)(xsi * ConcGridXRez) + 1;
                        int jk = (int)(eta * ConcGridYRez) + 1;
                        double[] depo_L = Program.Depo_conz[ik][jk];
                        double conc = masse * Pd1 * area_rez_fac;
                        lock (depo_L)
                        {
                            depo_L[SG_nteil] += conc;
                            depo_events++;
                        }
                        masse -= masse * Pd1;
                        depo_reflection_counter = -2; // block deposition 1 for 1 timestep

                        if (Program.ContinuousTraj)
                        {
                            if (masse < 1e-12) masse = 1e-12; // continuous mode: keep tiny residual mass
                        }
                        else
                        {
                            if (masse <= 0)
                            {
                                Program.Remove_MassZero++;//20260114
                                goto REMOVE_PARTICLE;
                            }
                        }
                    }
                }

                if (auszeit > 3600)
                {
                    distanceParticle += (float)Math.Sqrt(Program.Pow2(xcoord_nteil - xcoord_nteil_Prev) + Program.Pow2(ycoord_nteil - ycoord_nteil_Prev));
                    distanceParticle *= 0.5F;
                    if (distanceParticle < 0.02F)
                    {
                        goto REMOVE_PARTICLE;
                    }
                }

            MOVE_TO_CONCENTRATIONCALCULATION:

                //coordinate indices
                int iko = (int)(xsi * ConcGridXRez) + 1;
                int jko = (int)(eta * ConcGridYRez) + 1;

                if ((iko > Program.NXL) || (iko < 0) || (jko > Program.NYL) || (jko < 0))
                {
                    goto REMOVE_PARTICLE;
                }

                float zcoordRelative = 0;
                if (Program.BuildingsExist == true && Program.Topo == Consts.TerrainAvailable)
                {
                    //AHint contains the building height if terrain is available therefore the terrain height has to be corrected
                    zcoordRelative = zcoord_nteil - AHint + Program.CUTK[FFCellX][FFCellY];
                }
                else
                {
                    zcoordRelative = zcoord_nteil - AHint;
                }

                //Check for trapped particles
                {
                    int zko = (int)((zcoordRelative - AHint) * ConcGridZRez);
                    if (iko == ConcCellPrevX && jko == ConcCellPrevY && zko == ConcCellPrevZ)
                    {
                        ConcCellTime += idt;
                        if (ConcCellTime > ConcCellTimeMax - 10)
                        {
                            zcoord_nteil += 2;
                            if (ConcCellTime > ConcCellTimeMax)
                            {
                                goto REMOVE_PARTICLE;
                            }
                        }
                    }
                    else
                    {
                        ConcCellPrevX = iko;
                        ConcCellPrevY = jko;
                        ConcCellPrevZ = zko;
                        ConcCellTime = 0;
                    }
                }

                for (int II = 0; II < kko.Length; II++)
                {
                    float slice = (zcoordRelative - Program.HorSlices[II]) * ConcGridZRez;
                    if (Math.Abs(slice) > Int16.MaxValue)
                    {
                        goto REMOVE_PARTICLE;
                    }

                    kko[II] = (int) slice;
                }

                //decay rate
                if (decay_rate > 0)
                {
                    masse *= Math.Exp(-decay_rate * idt);
                }

                // compute Wet deposition
                if (Program.WetDeposition == true && Program.WetDepoRW != 0)
                {
                    double epsilonW_Masse = masse * Program.WetDepoRW * idt; // deposited mass
                    bool ex = false;

                    if (epsilonW_Masse > masse)
                    {
                        epsilonW_Masse = masse;
                        ex = true;
                    }
                    masse -= epsilonW_Masse;

                    double[] depo_L = Program.Depo_conz[iko][jko];
                    double conc = epsilonW_Masse * area_rez_fac;
                    lock (depo_L)
                    {
                        depo_L[SG_nteil] += conc;
                        depo_events++;//20251223
                    }

                    if (ex)
                    {
                        goto REMOVE_PARTICLE;
                    }
                }

                if (Deposition_type < Consts.DepoOnly) // compute concentrations for this particle
                {
                    //receptor concentrations
                    if ((Program.ReceptorsAvailable) && (reflexion_flag == Consts.ParticleNotReflected))
                    {
                        for (int irec = 1; irec < ReceptorConcentration.Length; irec++)
                        {
                            // if a receptor is inside or nearby a building, use the raster grid concentration inside or nearby the building
                            if (Program.ReceptorNearbyBuilding[irec])
                            {
                                if (iko == Program.ReceptorIInd[irec] &&
                                    jko == Program.ReceptorJInd[irec])
                                {
                                    float slice = (zcoord_nteil - AHint - Program.ReceptorZ[irec]) * ConcGridZRez;
                                    if ((int) slice == 0)
                                    {
                                        ReceptorConcentration[irec] += idt * masse;
                                    }
                                }
                            }
                            else // use the concentration at the receptor position x +- GralDx/2 and receptor y +- GralDy/2
                            {
                                if (irec == 1 && nteil == 1 && Program.IWET == 1)
                                {
                                    double dx = Math.Abs(xcoord_nteil - Program.ReceptorX[irec]);
                                    double dy = Math.Abs(ycoord_nteil - Program.ReceptorY[irec]);
                                    double dz = (zcoord_nteil - AHint - Program.ReceptorZ[irec]) * ConcGridZRez;
                                    System.Console.WriteLine($"[DBG] gate dx={dx:F2} dy={dy:F2} dz={dz:F2} refl={reflexion_flag}");
                                }
                                if (Math.Abs(xcoord_nteil - Program.ReceptorX[irec]) < ConcGridXHalf &&
                                    Math.Abs(ycoord_nteil - Program.ReceptorY[irec]) < ConcGridYHalf)
                                {
                                    float slice = (zcoord_nteil - AHint - Program.ReceptorZ[irec]) * ConcGridZRez;
                                    if ((int) slice == 0)
                                    {
                                        ReceptorConcentration[irec] += idt * masse;
                                    }
                                }
                            }
                        }
                    }

                    //compute 2D - concentrations
                    for (int II = 0; II < kko.Length; II++)
                    {
                        if ((kko[II] == 0) && (reflexion_flag == Consts.ParticleNotReflected))
                        {
                            float[] conz3d_L = Program.Conz3d[iko][jko][II];
                            double conc = masse * idt;
                            lock (conz3d_L)
                            {
                                conz3d_L[SG_nteil] += (float)conc;
                            }
                        }
                    }

                    //count particles in case of non-steady-state dispersion for the 3D concentration file
                    if (ISTATIONAER == Consts.TransientMode && Program.WriteVerticalConcentration)
                    {
                        if (reflexion_flag == Consts.ParticleNotReflected)
                        {
                            TransientGridZ = TransientConcentration.BinarySearchTransient(zcoord_nteil - AHint);
                            TransientGridX = (int)(xsi * FFGridXRez) + 1;
                            TransientGridY = (int)(eta * FFGridYRez) + 1;

                            float[] conzsum_L = Program.ConzSsum[TransientGridX][TransientGridY];
                            double conc = masse * Program.GridVolume * idt / (Area_cart * Program.DZK_Trans[TransientGridZ]);
                            lock (conzsum_L)
                            {
                                conzsum_L[TransientGridZ] += (float)conc;
                            }
                        }
                    }

                    //ODOUR Dispersion requires the concentration fields above and below the acutal layer
                    if (Program.Odour == true)
                    {
                        for (int II = 0; II < kko.Length; II++)
                        {
                            float slice = (zcoordRelative - (Program.HorSlices[II] + Program.GralDz)) * ConcGridZRez;
                            kko[II] = (int) (Math.Min(int.MaxValue, slice));
                        }
                        for (int II = 0; II < kko.Length; II++)
                        {
                            if ((kko[II] == 0) && (reflexion_flag == Consts.ParticleNotReflected))
                            {
                                float[] conz3dp_L = Program.Conz3dp[iko][jko][II];
                                double conc = masse * idt;
                                lock (conz3dp_L)
                                {
                                    conz3dp_L[SG_nteil] += (float)conc;
                                }
                            }
                        }
                        for (int II = 0; II < kko.Length; II++)
                        {
                            float slice = (zcoordRelative - (Program.HorSlices[II] - Program.GralDz)) * ConcGridZRez;
                            kko[II] = (int)(Math.Min(int.MaxValue, slice));
                        }
                        for (int II = 0; II < kko.Length; II++)
                        {
                            if ((kko[II] == 0) && (reflexion_flag == Consts.ParticleNotReflected))
                            {
                                float[] conz3dm_L = Program.Conz3dm[iko][jko][II];
                                double conc = masse * idt;
                                lock (conz3dm_L)
                                {
                                    conz3dm_L[SG_nteil] += (float)conc;
                                }
                            }
                        }
                    }
                } // compute concentrations
            } //END OF TIME-LOOP FOR PARTICLES

        if (!endOfStep)//20260114
        {
            Program.Remove_MaxLoops++; 
        }

        if (Program.ContinuousTraj && endOfStep)//20260113
        {
            Program.Keep_EndOfStep++;//20260114
            Program.Xcoord[nteil] = xcoord_nteil;
            Program.YCoord[nteil] = ycoord_nteil;
            Program.ZCoord[nteil] = zcoord_nteil;
            Program.ParticleMass[nteil] = masse;
            return;
        }

        REMOVE_PARTICLE:
        //20260129
        if (Program.ContinuousTraj && Program.Keep_Ground > 0 && Program.Remove_Other < 5)
        {
            System.Console.WriteLine("[DBG] removed after keep_ground");
        }

            if (Program.LogLevel > Consts.LogLevelOff) // additional log output
            {
                #region log_output
                if (reflexion_number > 400000)
                {
                    lock (Program.LogReflexions.SyncRoot)
                    {
                        Program.LogReflexions[5]++;
                    }
                }
                else if (reflexion_number > 100000)
                {
                    lock (Program.LogReflexions.SyncRoot)
                    {
                        Program.LogReflexions[4]++;
                    }
                }
                else if (reflexion_number > 50000)
                {
                    lock (Program.LogReflexions.SyncRoot)
                    {
                        Program.LogReflexions[3]++;
                    }
                }
                else if (reflexion_number > 10000)
                {
                    lock (Program.LogReflexions.SyncRoot)
                    {
                        Program.LogReflexions[2]++;
                    }
                }
                else if (reflexion_number > 5000)
                {
                    lock (Program.LogReflexions.SyncRoot)
                    {
                        Program.LogReflexions[1]++;
                    }
                }
                else if (reflexion_number > 500)
                {
                    lock (Program.LogReflexions.SyncRoot)
                    {
                        Program.LogReflexions[0]++;
                    }
                }

                if (timestep_number > 100000000)
                {
                    lock (Program.Log_Timesteps.SyncRoot)
                    {
                        Program.Log_Timesteps[5]++;
                    }
                }
                else if (timestep_number > 10000000)
                {
                    lock (Program.Log_Timesteps.SyncRoot)
                    {
                        Program.Log_Timesteps[4]++;
                    }
                }
                else if (timestep_number > 1000000)
                {
                    lock (Program.Log_Timesteps.SyncRoot)
                    {
                        Program.Log_Timesteps[3]++;
                    }
                }
                else if (timestep_number > 100000)
                {
                    lock (Program.Log_Timesteps.SyncRoot)
                    {
                        Program.Log_Timesteps[2]++;
                    }
                }
                else if (timestep_number > 10000)
                {
                    lock (Program.Log_Timesteps.SyncRoot)
                    {
                        Program.Log_Timesteps[1]++;
                    }
                }
                else if (timestep_number > 1000)
                {
                    lock (Program.Log_Timesteps.SyncRoot)
                    {
                        Program.Log_Timesteps[0]++;
                    }
                }

                #endregion log_output
            }

            // Add local receptor concentrations to receptor array and store maximum concentration part for each receptor
            if (Program.ReceptorsAvailable)
            {
                for (int ii = 1; ii < ReceptorConcentration.Length; ii++)
                {
                    double recconc = ReceptorConcentration[ii];
                    if (recconc > 0)
                    {
                        double[] recmit_L = Program.ReceptorConc[ii];
                        lock (recmit_L)
                        {
                            recmit_L[SG_nteil] += recconc;
                        }

                        if (recconc > Program.ReceptorParticleMaxConc[ii][SG_nteil])
                        {
                            Interlocked.Exchange(ref Program.ReceptorParticleMaxConc[ii][SG_nteil], recconc);
                        }
                    }
                }
            }
            Program.Remove_Other++;//20260114
            Program.ParticleSource[nteil] = 0;//20260113
            Program.ParticleMass[nteil] = 0;//20260113
            return;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
        static float Erf(float x) // Error function
        {
            // constants
            const float a1 = 0.254829592F;
            const float a2 = -0.284496736F;
            const float a3 = 1.421413741F;
            const float a4 = -1.453152027F;
            const float a5 = 1.061405429F;
            const float p = 0.3275911F;

            // Save the sign of x
            int sign = 1;
            if (x < 0)
            {
                sign = -1;
            }

            x = Math.Abs(x);

            // A&S formula 7.1.26
            float t = 1 / (1 + p * x);
            float y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * MathF.Exp(-x * x);
            return sign * y;
        }

        /// <summary>
        ///Calculate the deposition probability according to VDI 3945
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
        public static float Pd(float varw, float vsed, float vdep, int Deposition_Type, int IndexI, int IndexJ)
        {
            varw = MathF.Sqrt(varw);
            // inside a vegetation area?
            if (Program.COV[IndexI][IndexJ] > 0)
            {
                if (Deposition_Type == 2) // Deposition only - PM30 or larger
                {
                    vdep *= (1F + Program.VegetationDepoVelFactors.VelPMxxFact * Program.COV[IndexI][IndexJ]); // default: up to * 3 for large particles if COV[][] > 0
                }
                else
                {
                    vdep *= (1F + Program.VegetationDepoVelFactors.VelGasFact * Program.COV[IndexI][IndexJ]); // default: up to * 1.5 for small particles and gases if COV[][] > 0
                }
            }
            // compute deposition probability according to VDI 3945 for this particle
            float Ks = vsed / (sqrt2F * varw);
            float Fs = sqrtPiF * Ks + MathF.Exp(-MathF.Pow(Ks, 2)) / (1 + Erf(Ks));
            float Pd = 3 * sqrtPiF * 3600 / Program.TAUS * vdep / (sqrt2F * varw * Fs + sqrtPiF * vdep);
            if (Pd < 1.0)
            {
                return Pd;
            }
            else
            {
                return 1;
            }
        }

    }
}


