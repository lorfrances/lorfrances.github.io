const kB =  1.380649e-23;             // global constant
const hbar = 1.054571817e-34;
const m0 = 9.1093837e-31;
const e_c = 1.60217663e-19;
const eps0 = 8.8541878188e-12;
const Fieldy = 0.0;
const Fieldz = 0.0;
const T = 300;
const NumValley = 3;
const NumPhysValley = 5;
const NumPar = 5000;
const NumEnergyLevel = 1000;
const NumEnergyLevelStat = 500;
const EnergyStep = 1.0e-3; //eV
const EnergyStepStat = 2.0e-3;
const Dt = 5.0e-15;
const { mat } = initializeMat();
const scatt = ScatteringTable(mat);
const maxScatt = Math.max(...scatt.PolarGAbs) + Math.max(...scatt.PolarGEmit) + Math.max(...scatt.IVGXAbs) + Math.max(...scatt.IVGXEmit) + Math.max(...scatt.IVGLAbs) + Math.max(...scatt.IVGLEmit);

const Tau0 = 1/maxScatt;
const { state } = initializePar(NumPar, mat);
const LatConst = 5.65e-10;

// drawing constants and display stuff
const BZ_bound = (2*Math.PI)/LatConst;
const kMax = BZ_bound;
const theCanvas = document.getElementById("theCanvas");
const theContext = theCanvas.getContext("2d");
const CanvasHist = document.getElementById("CanvasHist");
const VelHistCanvas = document.getElementById("CanvasVel");
const HistContext = CanvasHist.getContext("2d")
const VelContext = VelHistCanvas.getContext("2d")
const BarWidth_v = (VelHistCanvas.width-20)/200
const BarWidth = (CanvasHist.width-20)/NumEnergyLevelStat
const BarScale = (CanvasHist.height-20)/(NumPar*0.05)
const kScale = theCanvas.width/(2*kMax);
const stepsFrame = 10;

const VelMax = 1.5e6 // m/s
const VelStep = 15000// m/s
const NumVelStep = 100


HistContext.strokeStyle = "black"; // line color
HistContext.lineWidth = 1;        // line thickness
HistContext.beginPath();
HistContext.moveTo(20, 100);  // start point
HistContext.lineTo(520, 100); // end point
HistContext.stroke();
HistContext.beginPath();
HistContext.moveTo(20, 100);  // start point
HistContext.lineTo(20, 10); // end point
HistContext.stroke();
HistContext.font= "16px Arial";
HistContext.fillText("Energy", 240, 117);
HistContext.fillText("N", 0, 50);

VelContext.strokeStyle = "black"; // line color
VelContext.lineWidth = 1;        // line thickness
VelContext.beginPath();
VelContext.moveTo(20, 100);  // start point
VelContext.lineTo(220, 100); // end point
VelContext.stroke();
VelContext.beginPath();
VelContext.moveTo(20, 100);  // start point
VelContext.lineTo(20, 10); // end point
VelContext.stroke();
VelContext.font= "16px Arial";
VelContext.fillText("Velocity", 85, 117);
VelContext.fillText("N", 0, 50);

// Standard Normal variate using Box-Muller transform.
function gaussianRandom(mean=0, stdev=1) {
    const u = 1 - Math.random(); // Converting [0,1) to (0,1]
    const v = Math.random();
    const z = Math.sqrt( -2.0 * Math.log( u ) ) * Math.cos( 2.0 * Math.PI * v );
    // Transform to the desired mean and standard deviation:
    return z * stdev + mean;
}

function initializeMat() {
    const mass = new Float64Array(NumValley); //density of states effective mass (for scattering)
    const Eg = new Float64Array(NumValley);

    Eg[0] = 1.42;
    Eg[1] = 1.71;
    Eg[2] = 1.90;

    // density of states masses
    mass[0] = 0.067;
    mass[1] = 0.22;
    mass[2] = 0.40;
    const m_X_l = 1.3;
    const m_X_t = 0.23;
    const m_L_l = 1.9;
    const m_L_t = 0.075;
    // anisotropic masses and transformation
    // Gamma valley
    const t_G = math.matrix([[1/Math.sqrt(mass[0]), 0.0, 0.0], [0.0, 1/Math.sqrt(mass[0]), 0.0], [0.0, 0.0, 1/Math.sqrt(mass[0])]])
    const t_G_inv = math.inv(t_G)
    const m_G_r = math.flatten(math.matrix([[1/mass[0], 0.0, 0.0], [0.0, 1/mass[0], 0.0], [0.0, 0.0, 1/mass[0]]])).toArray();

    // L valleys
    const R_L111 = math.matrix([[1/Math.sqrt(3), 1/Math.sqrt(3), 1/Math.sqrt(3)],[-1/Math.sqrt(2), 1/Math.sqrt(2), 0], [-1/Math.sqrt(6), -1/Math.sqrt(6), Math.sqrt(2)/Math.sqrt(3)]]);
    const r_L111 = math.flatten(R_L111).toArray(); 
    const s_L = math.matrix([[1/Math.sqrt(m_L_l), 0.0, 0.0], [0.0, 1/Math.sqrt(m_L_t), 0.0], [0.0, 0.0, 1/Math.sqrt(m_L_t)]]);
    const t_L111 = math.multiply(s_L, R_L111);
    const t_L111_inv = math.inv(t_L111)
    const m_L111_inv = math.matrix([[1/m_L_l, 0.0, 0.0], [0.0, 1/m_L_t, 0.0], [0.0, 0.0, 1/m_L_t]])
    const m_L111_r = math.flatten(math.multiply(math.transpose(R_L111),m_L111_inv,R_L111)).toArray();

    const R_L1_11 = math.matrix([[-1/Math.sqrt(3), 1/Math.sqrt(3), 1/Math.sqrt(3)],[-1/Math.sqrt(2), -1/Math.sqrt(2), 0], [1/Math.sqrt(6), -1/Math.sqrt(6), Math.sqrt(2)/Math.sqrt(3)]]);
    const r_L1_11 = math.flatten(R_L1_11).toArray(); 
    const t_L1_11 = math.multiply(s_L, R_L1_11);
    const t_L1_11_inv = math.inv(t_L1_11);
    const m_L1_11_r = math.flatten(math.multiply(math.transpose(R_L1_11),m_L111_inv,R_L1_11)).toArray();

    const R_L11_1 = math.matrix([[1/Math.sqrt(3), -1/Math.sqrt(3), 1/Math.sqrt(3)],[1/Math.sqrt(2), 1/Math.sqrt(2), 0], [-1/Math.sqrt(6), 1/Math.sqrt(6), Math.sqrt(2)/Math.sqrt(3)]]);
    const r_L11_1 = math.flatten(R_L11_1).toArray(); 
    const t_L11_1 = math.multiply(s_L, R_L11_1);
    const t_L11_1_inv = math.inv(t_L11_1);
    const m_L11_1_r = math.flatten(math.multiply(math.transpose(R_L11_1),m_L111_inv,R_L11_1)).toArray()

    const R_L111_ = math.matrix([[1/Math.sqrt(3), 1/Math.sqrt(3), -1/Math.sqrt(3)],[-1/Math.sqrt(2), 1/Math.sqrt(2), 0], [1/Math.sqrt(6), 1/Math.sqrt(6), Math.sqrt(2)/Math.sqrt(3)]]);
    const r_L111_ = math.flatten(R_L111_).toArray(); 
    const t_L111_ = math.multiply(s_L, R_L111_);
    const t_L111__inv = math.inv(t_L111_);
    const m_L111__r = math.flatten(math.multiply(math.transpose(R_L111_),m_L111_inv,R_L111_)).toArray()


    // X valleys
    const t_X100 = math.matrix([[1/Math.sqrt(m_X_l), 0.0, 0.0], [0.0, 1/Math.sqrt(m_X_t), 0.0], [0.0, 0.0, 1/Math.sqrt(m_X_t)]]);
    const t_X010 = math.matrix([[1/Math.sqrt(m_X_t), 0.0, 0.0], [0.0, 1/Math.sqrt(m_X_l), 0.0], [0.0, 0.0, 1/Math.sqrt(m_X_t)]]);
    const t_X001 = math.matrix([[1/Math.sqrt(m_X_t), 0.0, 0.0], [0.0, 1/Math.sqrt(m_X_t), 0.0], [0.0, 0.0, 1/Math.sqrt(m_X_l)]]);
    const m_X100_r = math.flatten(math.matrix([[1/(m_X_l), 0.0, 0.0], [0.0, 1/(m_X_t), 0.0], [0.0, 0.0, 1/(m_X_t)]])).toArray();
    const m_X010_r = math.flatten(math.matrix([[1/(m_X_t), 0.0, 0.0], [0.0, 1/(m_X_l), 0.0], [0.0, 0.0, 1/(m_X_t)]])).toArray();
    const m_X001_r = math.flatten(math.matrix([[1/(m_X_t), 0.0, 0.0], [0.0, 1/(m_X_t), 0.0], [0.0, 0.0, 1/(m_X_l)]])).toArray();

    const t_X100_inv = math.inv(t_X100)
    const t_X010_inv = math.inv(t_X010)
    const t_X001_inv = math.inv(t_X001)
    
    const T_G = math.flatten(t_G).toArray();
    const T_X100 = math.flatten(t_X100).toArray();
    const T_X010 = math.flatten(t_X010).toArray();
    const T_X001 = math.flatten(t_X001).toArray();
    const T_L111 = math.flatten(t_L111).toArray();
    const T_L1_11 = math.flatten(t_L1_11).toArray();
    const T_L1_11_inv = math.flatten(t_L1_11_inv).toArray();
    const T_L11_1 = math.flatten(t_L11_1).toArray();
    const T_L11_1_inv = math.flatten(t_L11_1_inv).toArray();
    const T_L111_ = math.flatten(t_L111_).toArray();
    const T_L111__inv = math.flatten(t_L111__inv).toArray();
    const T_G_inv = math.flatten(t_G_inv).toArray();
    const T_X100_inv = math.flatten(t_X100_inv).toArray();
    const T_X010_inv = math.flatten(t_X010_inv).toArray();
    const T_X001_inv = math.flatten(t_X001_inv).toArray();
    const T_L111_inv = math.flatten(t_L111_inv).toArray();  

    const eps_s = 12.90;
    const eps_inf = 10.89;
    const EpO = 0.036;
    const dens = 5320;
    return {
        mat: {mass, eps_s, eps_inf, EpO, dens, Eg, T_X100, T_X010, T_X001, T_G, T_L111, T_G_inv, T_X100_inv, T_X010_inv, T_X001_inv, T_L111_inv, T_L1_11, T_L1_11_inv, T_L11_1, T_L11_1_inv, T_L111_, T_L111__inv, r_L111,  r_L1_11, r_L11_1, r_L111_, m_G_r, m_L111_r, m_L1_11_r, m_L11_1_r, m_L111__r, m_X100_r, m_X010_r, m_X001_r}
    }
}

function initializePar(NumPar, mat) {   
    const kx  = new Float64Array(NumPar);    // initialize electron vectors
    const ky  = new Float64Array(NumPar);
    const kz  = new Float64Array(NumPar);
    const parTau = new Float64Array(NumPar);
    const valley = new Float64Array(NumPar);
    const  physValley = new Float64Array(NumPar);
    const sigma_k = (mat.mass[0]*m0/hbar)*Math.sqrt(kB*T/(mat.mass[0]*m0))
    for (let i = 0; i < NumPar; i++) {
        kx[i] = gaussianRandom(0, sigma_k);
        ky[i] = gaussianRandom(0, sigma_k);
        kz[i] = gaussianRandom(0, sigma_k);
        parTau[i] = -Tau0*Math.log(Math.random());
        valley[i] = 0;
        physValley[i] = 0;
    }
    return {
           state: {kx, ky, kz, parTau, valley, physValley},
           
    }
}

// function acousticRate(E_k, mat) {
//     const DefPot = 7.0;
//     const W_a = (((2*mat.mass[0]*m0)**1.5 * kB*T * (e_c*DefPot)*(e_c*DefPot))/(2*Math.PI*5320*5000*5000*hbar**4))*Math.sqrt(E_k)
//     return W_a
// }

function polarRate(E_k, Abs, mat) {
    const Omega_pO = mat.EpO*e_c/hbar;
    let E_f;
    if (Abs == 1) {
        E_f = E_k + mat.EpO*e_c;
    }
    else {
        E_f = E_k - mat.EpO*e_c;
    }
    if (E_f/e_c > 1.0e-6) {
        const gamma = E_k;
        const gamma_f = E_f;
        const A = 4.0;
        const B = 0.0;
        const C = 4.0;
        const F0 = 1/C *(A * Math.log(Math.abs((Math.sqrt(gamma) + Math.sqrt(gamma_f))/(Math.sqrt(gamma)-Math.sqrt(gamma_f)))) + B);
        const Const = e_c*e_c * Math.sqrt(mat.mass[0]*m0) * Omega_pO/(4*Math.PI*Math.sqrt(2)*hbar) * ((1/(eps0*mat.eps_inf))-(1/(eps0*mat.eps_s)));
        const f_p = 1/(Math.exp(mat.EpO*e_c/(kB*T))-1)
        let W;
        if (Abs == 1) {
            W = Const*(1/(Math.sqrt(gamma)))*F0*f_p;
        }
        else {
            W = Const*(1/(Math.sqrt(gamma)))*F0*(f_p + 1);
        }
        return W;
    }
    else {
        return 0.0;
    }

}

function intervalleyRate(E_k, Z_f, mat, valley_i, valley_f, Ep, DefPot, Abs) {
    const Omega_piv = Ep*e_c/hbar;
    const f_p = 1/(Math.exp((Ep*e_c)/(kB*T)) - 1);
    const Valley_sep = mat.Eg[valley_f] - mat.Eg[valley_i];
    let E_f;
    if (Abs === 1) {
        E_f = E_k + e_c*Ep - Valley_sep*e_c;
    }
    else {
        E_f = E_k - e_c*Ep - Valley_sep*e_c;
    }
    const Const =    Z_f * Math.sqrt(2) * (m0*mat.mass[valley_f])**1.5 * (e_c*DefPot)**2 / (Math.PI * mat.dens * Omega_piv * hbar**3); 
    let W_IV;
    if (E_f/e_c > 1e-6) {
        if (Abs === 1) {
            W_IV = Const*E_f**0.5*f_p
        }
        else {
            W_IV = Const*E_f**0.5*(f_p+1)
        }
        return W_IV;
    }
    else {
        return W_IV = 0.0;
    }


}

function ScatteringTable(mat) {
    const PolarGAbs = new Float64Array(NumEnergyLevel);
    const PolarGEmit = new Float64Array(NumEnergyLevel);
    const IVGLAbs = new Float64Array(NumEnergyLevel);
    const IVGLEmit = new Float64Array(NumEnergyLevel);
    const IVGXAbs = new Float64Array(NumEnergyLevel);
    const IVGXEmit = new Float64Array(NumEnergyLevel);

    const IVLGAbs = new Float64Array(NumEnergyLevel);
    const IVLGEmit = new Float64Array(NumEnergyLevel);
    const IVLXAbs = new Float64Array(NumEnergyLevel);
    const IVLXEmit = new Float64Array(NumEnergyLevel);
    const IVLLAbs = new Float64Array(NumEnergyLevel);
    const IVLLEmit = new Float64Array(NumEnergyLevel);

    const IVXGAbs = new Float64Array(NumEnergyLevel);
    const IVXGEmit = new Float64Array(NumEnergyLevel);
    const IVXLAbs = new Float64Array(NumEnergyLevel);
    const IVXLEmit = new Float64Array(NumEnergyLevel);
    const IVXXAbs = new Float64Array(NumEnergyLevel);
    const IVXXEmit = new Float64Array(NumEnergyLevel);

    for (let j = 0; j < NumEnergyLevel; j++) {
        const E_k = (j+1)*EnergyStep*e_c
        PolarGAbs[j] = polarRate(E_k,1,mat);
        PolarGEmit[j] = polarRate(E_k,0,mat);
        IVGLAbs[j] = intervalleyRate(E_k, 4, mat, 0, 1, 0.028, 1e11, 1)
        IVGLEmit[j] = intervalleyRate(E_k, 4, mat, 0, 1, 0.028, 1e11, 0)
        IVGXAbs[j] = intervalleyRate(E_k, 3, mat, 0, 2, 0.030, 1e11, 1)
        IVGXEmit[j] = intervalleyRate(E_k, 3, mat, 0, 2, 0.030, 1e11, 0)

        IVLGAbs[j] = intervalleyRate(E_k, 1, mat, 1, 0, 0.028, 1e11, 1)
        IVLGEmit[j] = intervalleyRate(E_k, 1, mat, 1, 0, 0.028, 1e11, 0)
        IVLXAbs[j] = intervalleyRate(E_k, 3, mat, 1, 2, 0.030, 1e11, 1)
        IVLXEmit[j] = intervalleyRate(E_k, 3, mat, 1, 2, 0.030, 1e11, 0)
        IVLLAbs[j] = intervalleyRate(E_k, 3, mat, 1, 1, 0.030, 0.7e11, 1)
        IVLLEmit[j] = intervalleyRate(E_k, 3, mat, 1, 1, 0.030, 0.7e11, 0)

        IVXGAbs[j] = intervalleyRate(E_k, 1, mat, 2, 0, 0.030, 1e11, 1)
        IVXGEmit[j] = intervalleyRate(E_k, 1, mat, 2, 0, 0.030, 1e11, 0)
        IVXLAbs[j] = intervalleyRate(E_k, 4, mat, 2, 1, 0.030, 1e11, 1)
        IVXLEmit[j] = intervalleyRate(E_k, 4, mat, 2, 1, 0.030, 1e11, 0)
        IVXXAbs[j] = intervalleyRate(E_k, 2, mat, 2, 2, 0.030, 0.7e11, 1)
        IVXXEmit[j] = intervalleyRate(E_k, 2, mat, 2, 2, 0.030, 0.7e11, 0)

    }
    return {PolarGAbs, PolarGEmit, IVGLAbs, IVGLEmit, IVGXAbs, IVGXEmit, IVLGAbs, IVLGEmit, IVLXAbs, IVLXEmit, IVLLAbs, IVLLEmit, IVXGAbs, IVXGEmit, IVXLAbs, IVXLEmit, IVXXAbs, IVXXEmit};
}

function LoopParticles(state, scatt, Fieldx) {
    for (let i = 0; i < NumPar; i++) {
        let remainTime = Dt
        let loop = 1
        let carryTime = 0.0
        let newTau
        let driftTime
        while (remainTime > 0.0) {
            if (loop === 1) {
                newTau = state.parTau[i]
            }
            else {
                newTau = -Tau0*Math.log(Math.random())
            }

            if (newTau <= remainTime) {
                driftTime = newTau
                remainTime = remainTime - newTau
                state.kx[i] = state.kx[i] + driftTime*(-e_c*Fieldx/hbar)
                state.ky[i] = state.ky[i] + driftTime*(-e_c*Fieldy/hbar)
                state.kz[i] = state.kz[i] + driftTime*(-e_c*Fieldz/hbar)
                state = scatter(state, mat, i, scatt)
            }
            else {
                driftTime = remainTime
                carryTime = newTau - remainTime // carry over to next time step
                remainTime = 0.0
                state.kx[i] = state.kx[i] + driftTime*(-e_c*Fieldx/hbar)
                state.ky[i] = state.ky[i] + driftTime*(-e_c*Fieldy/hbar)
                state.kz[i] = state.kz[i] + driftTime*(-e_c*Fieldz/hbar)               
            }
            loop = loop + 1
        }
        state.parTau[i] = carryTime
    }
    return state
}

function scatter(state, mat, par, scatt) {
    let T_ij;
    if (state.physValley[par] === 0) {
        T_ij = mat.T_G;
    }
    else if (state.physValley[par] === 1) {
        T_ij = mat.T_L111;
    }
    else if (state.physValley[par] === 2){
        T_ij = mat.T_L1_11;
    }
    else if (state.physValley[par] === 3){
        T_ij = mat.T_L11_1;
    }
    else if (state.physValley[par] === 4){
        T_ij = mat.T_L111_;
    }   
    else if (state.physValley[par] === 5) {
        T_ij = mat.T_X100;
    }
    else if (state.physValley[par] === 6) {
        T_ij = mat.T_X010;
    }
    else if (state.physValley[par] === 7) {
        T_ij = mat.T_X001;
    }

    const w_x = state.kx[par] * T_ij[0] + state.ky[par] * T_ij[1] + state.kz[par] * T_ij[2];
    const w_y = state.kx[par] * T_ij[3] + state.ky[par] * T_ij[4] + state.kz[par] * T_ij[5];
    const w_z = state.kx[par] * T_ij[6] + state.ky[par] * T_ij[7] + state.kz[par] * T_ij[8];
    const ParEn = (hbar*hbar * (w_x*w_x + w_y*w_y + w_z*w_z))/(2*m0*e_c) ;// eV energy
    let EnInx = Math.round(ParEn/EnergyStep);
    if (EnInx === 0) {
        EnInx = 1;
    }
    else if (EnInx >= NumEnergyLevel) {
        EnInx = NumEnergyLevel
    }
    let nScattMech;
    let atEMechs;
    if (state.valley[par] === 0) {
        atEMechs = [
            scatt.PolarGAbs[EnInx-1],
            scatt.PolarGEmit[EnInx-1],
            scatt.IVGLAbs[EnInx-1],
            scatt.IVGLEmit[EnInx-1],
            scatt.IVGXAbs[EnInx-1],
            scatt.IVGXEmit[EnInx-1],
            ];
        nScattMech = 6
    }
    else if (state.valley[par] === 1) {
        atEMechs = [
            scatt.IVLGAbs[EnInx-1],
            scatt.IVLGEmit[EnInx-1],
            scatt.IVLXAbs[EnInx-1],
            scatt.IVLXEmit[EnInx-1],
            scatt.IVLLAbs[EnInx-1],
            scatt.IVLLEmit[EnInx-1],
            ];
       nScattMech = 6
    }
    else if (state.valley[par] === 2) {
        atEMechs = [
            scatt.IVXGAbs[EnInx-1],
            scatt.IVXGEmit[EnInx-1],
            scatt.IVXLAbs[EnInx-1],
            scatt.IVXLEmit[EnInx-1],
            scatt.IVXXAbs[EnInx-1],
            scatt.IVXXEmit[EnInx-1],
            ];
       nScattMech = 6
    }
    const GammaInx = Math.random()*(1/Tau0);
    let SelectMech = -1
    let RateLo = 0.0
    let RateHi = 0.0
    for (let j = 0; j < nScattMech; j++) {
        RateLo = RateHi;
        RateHi = atEMechs[j] + RateHi;
        if ((GammaInx >= RateLo) && (GammaInx < RateHi)) {
            SelectMech = j
            break
        }
    }
    if (state.valley[par] === 0) {
        if (SelectMech === 0) {
            state = updatePolar(par, ParEn, state, 1, mat);
        }
        else if (SelectMech === 1) {
            state = updatePolar(par, ParEn, state, 0, mat);
        }
        else if (SelectMech === 2) {
            state = updateIV(par, ParEn, state, 0.028, 1, 1, mat, 4);
        }
        else if (SelectMech === 3) {
            state = updateIV(par, ParEn, state, 0.028, 1, 0, mat, 4);
        }
        else if (SelectMech === 4) {
            state = updateIV(par, ParEn, state, 0.030, 2, 1, mat, 3);
        }
        else if (SelectMech === 5) {
            state = updateIV(par, ParEn, state, 0.030, 2, 0, mat, 3);
        }
    }
    else if (state.valley[par] === 1) {
        if (SelectMech === 0) {
            state = updateIV(par, ParEn, state, 0.030, 0, 1, mat, 1);
        }
        else if (SelectMech === 1) {
            state = updateIV(par, ParEn, state, 0.030, 0, 0, mat, 1);
        }
        else if (SelectMech === 2) {
            state = updateIV(par, ParEn, state, 0.030, 2, 1, mat, 3);
        }
        else if (SelectMech === 3) {
            state = updateIV(par, ParEn, state, 0.030, 2, 0, mat, 3);
        }
        else if (SelectMech === 4) {
            state = updateIV(par, ParEn, state, 0.030, 1, 1, mat, 3);
        }
        else if (SelectMech === 5) {
            state = updateIV(par, ParEn, state, 0.030, 1, 0, mat, 3);
        }
    }
    else if (state.valley[par] === 2) {
        if (SelectMech === 0) {
            state = updateIV(par, ParEn, state, 0.030, 0, 1, mat, 1);
        }
        else if (SelectMech === 1) {
            state = updateIV(par, ParEn, state, 0.030, 0, 0, mat, 1);
        }
        else if (SelectMech === 2) {
            state = updateIV(par, ParEn, state, 0.030, 1, 1, mat, 4);
        }
        else if (SelectMech === 3) {
            state = updateIV(par, ParEn, state, 0.030, 1, 0, mat, 4);
        }
        else if (SelectMech === 4) {
            state = updateIV(par, ParEn, state, 0.030, 2, 1, mat, 2 );
        }
        else if (SelectMech === 5) {
            state = updateIV(par, ParEn, state, 0.030, 2, 0, mat, 2);
        }
    }
    return state;
}

function updateAcou(par, state) {
    const Wvk_p = Math.sqrt(state.kx[par]*state.kx[par] + state.ky[par]*state.ky[par] + state.kz[par]*state.kz[par]);
    const fi = 2*Math.PI*Math.random();
    const costh = 1-2*Math.random();
    const sinth = Math.sqrt(1 - costh*costh);
    state.kx[par] = Wvk_p*sinth*Math.cos(fi);
    state.ky[par] = Wvk_p*sinth*Math.sin(fi);
    state.kz[par] = Wvk_p*costh;
    return state;
}

function updateIV(par, ParEn, state, Ep, valley_f, Abs, mat, Z_f) {
    // handle anisotropy
    const Valley_sep = mat.Eg[valley_f] - mat.Eg[state.valley[par]];
    let NewEn;
    if (Abs === 1) {
        NewEn = (ParEn + Ep -  Valley_sep)*e_c;
    }    
    else {
        NewEn = (ParEn - Ep -  Valley_sep)*e_c;
    }

    let T_ij;
    if (valley_f === 0) {
        T_ij = mat.T_G_inv;
        state.physValley[par] = 0;
    }
    else if (valley_f === 1) {
        draw_valley = Math.floor(Math.random()*4) + 1;
        while (draw_valley === state.physValley[par]) {
            draw_valley = Math.floor(Math.random()*4) + 1;
        }
        state.physValley[par] = draw_valley;
        if (draw_valley === 1) {
            // pick 111
            T_ij = mat.T_L111_inv;
        }
        else if (draw_valley === 2) {
            // pick 1-11
            T_ij = mat.T_L1_11_inv;
        }
        else if (draw_valley === 3) {
            // pick 1-11
            T_ij = mat.T_L11_1_inv;
        }
        else if (draw_valley === 4) {
            // pick 1-11
            T_ij = mat.T_L111__inv;
        }
    }
    else if (valley_f === 2) { // if we end up in X valley, pick a subvalley rnadomly and uniformly bc we used dos effective mass to scatter
        draw_valley = Math.floor(Math.random()*3) + 5;
        while (draw_valley === state.physValley[par]) {
            draw_valley = Math.floor(Math.random()*3) + 5;
        }
        state.physValley[par] = draw_valley;
        if (draw_valley === 5) {
            // pick 100
            T_ij = mat.T_X100_inv;
        }
        else if (draw_valley === 6) {
            // pick 010
            T_ij = mat.T_X010_inv;
        }
        else {
            T_ij = mat.T_X001_inv;
        }
    }
    
    const w_mag = Math.sqrt(((2 * m0 * NewEn))/(hbar*hbar)) // 
    const fi = 2*Math.PI*Math.random();
    const costh = 1-2*Math.random();
    const sinth = Math.sqrt(1 - costh*costh);
    const w_x = w_mag*sinth*Math.cos(fi);
    const w_y = w_mag*sinth*Math.sin(fi);
    const w_z = w_mag*costh;
    state.kx[par] = w_x * T_ij[0] + w_y * T_ij[1] + w_z * T_ij[2];
    state.ky[par] = w_x * T_ij[3] + w_y * T_ij[4] + w_z * T_ij[5];
    state.kz[par] = w_x * T_ij[6] + w_y * T_ij[7] + w_z * T_ij[8];
    state.valley[par] = valley_f
    return state
}

function updatePolar(par, ParEn, state, Abs, mat) {
    let NewEn;
    if (Abs === 1){
        NewEn = (ParEn + mat.EpO)*e_c;
    }
    else {
        NewEn = (ParEn - mat.EpO)*e_c;
    }
    const Wvk_xy = Math.sqrt(state.kx[par]*state.kx[par] + state.ky[par]*state.ky[par]);
    const Wvk_xy_z = Math.sqrt(Wvk_xy*Wvk_xy + state.kz[par]*state.kz[par]);
    const costh0 = state.kz[par]/Wvk_xy_z;
    const sinth0 = Wvk_xy/Wvk_xy_z;
    const cosph0 = state.kx[par]/Wvk_xy;
    const sinph0 = state.ky[par]/Wvk_xy;   

    const Wvk_p = Math.sqrt(((2 * (mat.mass[0]*m0) * NewEn))/(hbar*hbar));
    const zeta = 2*Math.sqrt(ParEn*e_c*NewEn)/(ParEn*e_c+NewEn-2*Math.sqrt(ParEn*e_c*NewEn));
    const costh = ((zeta+1)-(2*zeta+1)**Math.random())/zeta;
    const sinth = Math.sqrt(1-costh*costh);
    const fai = 2*Math.PI*Math.random();
    const cosfai = Math.cos(fai);
    const sinfai = Math.sin(fai);
    const kxp = Wvk_p*sinth*cosfai;
    const kyp = Wvk_p*sinth*sinfai;
    const kzp = Wvk_p*costh;
    state.kx[par] = kxp*cosph0*costh0-kyp*sinph0+kzp*cosph0*sinth0;
    state.ky[par] = kxp*sinph0*costh0+kyp*cosph0+kzp*sinph0*sinth0;
    state.kz[par] = -kxp*sinth0+kzp*costh0;
    return state;
}

function statistics(state, totals) {
    for (let i = 0; i < NumPar; i++) {
        let T_ij;
        let M_ij;
        if (state.physValley[i] === 0) {
            T_ij = mat.T_G;
            M_ij = mat.m_G_r;
        }
        else if (state.physValley[i] === 1) {
            T_ij = mat.T_L111;
            M_ij = mat.m_L111_r;
        }
        else if (state.physValley[i] === 2){
            T_ij = mat.T_L1_11;
            M_ij = mat.m_L1_11_r;
        }
        else if (state.physValley[i] === 3){
            T_ij = mat.T_L11_1;
            M_ij = mat.m_L11_1_r;
        }
        else if (state.physValley[i] === 4){
            T_ij = mat.T_L111_;
            M_ij = mat.m_L111__r;
        }   
        else if (state.physValley[i] === 5) {
            T_ij = mat.T_X100;
            M_ij = mat.m_X100_r;
        }
        else if (state.physValley[i] === 6) {
            T_ij = mat.T_X010;
            M_ij = mat.m_X010_r
        }
        else if (state.physValley[i] === 7) {
            T_ij = mat.T_X001;
            M_ij = mat.m_X001_r
        }
        const w_x = state.kx[i] * T_ij[0] + state.ky[i] * T_ij[1] + state.kz[i] * T_ij[2];
        const w_y = state.kx[i] * T_ij[3] + state.ky[i] * T_ij[4] + state.kz[i] * T_ij[5];
        const w_z = state.kx[i] * T_ij[6] + state.ky[i] * T_ij[7] + state.kz[i] * T_ij[8];
        const v_x = (hbar/m0) * (M_ij[0]*state.kx[i] + M_ij[1]*state.ky[i] + M_ij[2]*state.kz[i]);
        const ParEn = (hbar*hbar * (w_x*w_x + w_y*w_y + w_z*w_z))/(2*m0*e_c) ;// eV energy
        const ParEnV = ParEn + (mat.Eg[state.valley[i]] - mat.Eg[0]);
        totals.En = totals.En + ParEn;
        totals.Vel = totals.Vel + v_x;
        let EnInx = Math.round(ParEnV/EnergyStepStat);
        if (EnInx === 0) {
            EnInx = 1;
        }
        else if (EnInx >= NumEnergyLevelStat){
            EnInx = NumEnergyLevelStat;
        }
        totals.EnHist[EnInx-1] = totals.EnHist[EnInx-1] + 1;

        let VelInx = Math.round(v_x/VelStep) + NumVelStep;
        if (VelInx <= 0) {
            VelInx = 0;
        }
        else if (VelInx >= 200){
            VelInx = 200;
        }
        totals.VelHist[VelInx] = totals.VelHist[VelInx] + 1;       
    }
    return totals;
}

function averages(totals) {
    const aveEnHist = new Float64Array(NumEnergyLevelStat);
    aveEnHist.fill(0)
    const aveVelHist = new Float64Array(200);
    aveVelHist.fill(0)
    const aveEn = totals.En/stepsFrame/NumPar;
    const aveVel = totals.Vel*1e-6*100/stepsFrame/NumPar;
    for (let j = 0; j < NumEnergyLevelStat; j++) {
        aveEnHist[j] = totals.EnHist[j]/stepsFrame;
    }
    for (let j = 0; j < 200; j++) {
        aveVelHist[j] = totals.VelHist[j]/stepsFrame;
    }
    return {aveEnHist, aveVelHist, aveEn, aveVel}
}

function drawVelHist(avgs){
    VelContext.clearRect(20,0,200,100)
    VelContext.fillStyle = "steelblue"
    for (let j = 0; j < 200; j++) {
        const x = 20 + j*BarWidth_v
        const y = (VelHistCanvas.height-20) - avgs.aveVelHist[j]*BarScale;
        VelContext.fillRect(x,y,BarWidth_v,avgs.aveVelHist[j]*BarScale);
    }
    veltoPrint = Math.round(avgs.aveVel*1000)/1000;
    VelContext.fillStyle = "black";
    VelContext.font= "20px Arial";
    VelContext.fillText(`<v> = ${veltoPrint.toFixed(1)}`, 25, 25)
    VelContext.font= "20px Arial";
    VelContext.fillText("x10 cm/s", 130, 25)
    VelContext.font= "13px Arial";
    VelContext.fillText("6", 160, 12)
}

function drawEnHist(avgs){
    HistContext.clearRect(20,0,500,100)
    HistContext.fillStyle = "steelblue"
    for (let j = 0; j < NumEnergyLevelStat; j++) {
        const x = 20 + j*BarWidth
        const y = (CanvasHist.height-20)- avgs.aveEnHist[j]*BarScale;
        HistContext.fillRect(x,y,BarWidth,avgs.aveEnHist[j]*BarScale);
    }
    const enetoPrint = Math.round(avgs.aveEn*1000)/1000;
    HistContext.fillStyle = "black";
    HistContext.font= "20px Arial";
    HistContext.fillText(`<E> = ${enetoPrint.toFixed(3)} eV`, 350, 25);

}

function drawKappa(state,mat) {
    
    theContext.clearRect(0,0,500,500)
    theContext.strokeStyle = "white"; // line color
    theContext.lineWidth = 1;        // line thickness
    theContext.beginPath();
    theContext.moveTo(0, 500/2);  // start point
    theContext.lineTo(500, 500/2); // end point
    theContext.stroke();
    
    theContext.beginPath();
    theContext.moveTo(500/2, 0);  // start point
    theContext.lineTo(500/2, 500); // end point
    theContext.stroke();    
    for (let i = 0; i < NumPar; i++) {
        const v = state.physValley[i];
        x = state.kx[i];
        y = state.ky[i];
        z = state.kz[i];
        if (v===1) {
            const l_max = Math.sqrt(3)*Math.PI/LatConst;
            x = state.kx[i] + Math.PI/LatConst;
            y = state.ky[i] + Math.PI/LatConst;
            z = state.kz[i] + Math.PI/LatConst;

            let l = x * mat.r_L111[0] + y * mat.r_L111[1] + z * mat.r_L111[2];
            let t1 = x * mat.r_L111[3] + y * mat.r_L111[4] + z * mat.r_L111[5];
            let t2 = x * mat.r_L111[6] + y * mat.r_L111[7] + z * mat.r_L111[8];
            let l_new;
            if (l >= l_max) {
                l_new = -l_max + (l - l_max);
                x = l_new*mat.r_L111[0] + t1*mat.r_L111[3] + t2*mat.r_L111[6];
                y = l_new*mat.r_L111[1] + t1*mat.r_L111[4] + t2*mat.r_L111[7];
                z = l_new*mat.r_L111[2] + t1*mat.r_L111[5] + t2*mat.r_L111[8];
            }
        }
        else if (v===2) {
            const l_max = Math.sqrt(3)*Math.PI/LatConst;
            x = state.kx[i] - Math.PI/LatConst;
            y = state.ky[i] + Math.PI/LatConst;
            z = state.kz[i] + Math.PI/LatConst;

            let l = x * mat.r_L1_11[0] + y * mat.r_L1_11[1] + z * mat.r_L1_11[2];
            let t1 = x * mat.r_L1_11[3] + y * mat.r_L1_11[4] + z * mat.r_L1_11[5];
            let t2 = x * mat.r_L1_11[6] + y * mat.r_L1_11[7] + z * mat.r_L1_11[8];
            let l_new;
            if (l >= l_max) {
                l_new = -l_max + (l - l_max);
                x = l_new*mat.r_L1_11[0] + t1*mat.r_L1_11[3] + t2*mat.r_L1_11[6];
                y = l_new*mat.r_L1_11[1] + t1*mat.r_L1_11[4] + t2*mat.r_L1_11[7];
                z = l_new*mat.r_L1_11[2] + t1*mat.r_L1_11[5] + t2*mat.r_L1_11[8];
            }
        }
        else if (v===3) {
            const l_max = Math.sqrt(3)*Math.PI/LatConst;
            x = state.kx[i] + Math.PI/LatConst;
            y = state.ky[i] - Math.PI/LatConst;
            z = state.kz[i] + Math.PI/LatConst;

            let l = x * mat.r_L11_1[0] + y * mat.r_L11_1[1] + z * mat.r_L11_1[2];
            let t1 = x * mat.r_L11_1[3] + y * mat.r_L11_1[4] + z * mat.r_L11_1[5];
            let t2 = x * mat.r_L11_1[6] + y * mat.r_L11_1[7] + z * mat.r_L11_1[8];
            let l_new;
            if (l >= l_max) {
                l_new = -l_max + (l - l_max);
                x = l_new*mat.r_L11_1[0] + t1*mat.r_L11_1[3] + t2*mat.r_L11_1[6];
                y = l_new*mat.r_L11_1[1] + t1*mat.r_L11_1[4] + t2*mat.r_L11_1[7];
                z = l_new*mat.r_L11_1[2] + t1*mat.r_L11_1[5] + t2*mat.r_L11_1[8];
            }
        }
        else if (v===4) {
            const l_max = Math.sqrt(3)*Math.PI/LatConst;
            x = state.kx[i] + Math.PI/LatConst;
            y = state.ky[i] + Math.PI/LatConst;
            z = state.kz[i] - Math.PI/LatConst;

            let l = x * mat.r_L111_[0] + y * mat.r_L111_[1] + z * mat.r_L111_[2];
            let t1 = x * mat.r_L111_[3] + y * mat.r_L111_[4] + z * mat.r_L111_[5];
            let t2 = x * mat.r_L111_[6] + y * mat.r_L111_[7] + z * mat.r_L111_[8];
            let l_new;
            if (l >= l_max) {
                l_new = -l_max + (l - l_max);
                x = l_new*mat.r_L111_[0] + t1*mat.r_L111_[3] + t2*mat.r_L111_[6];
                y = l_new*mat.r_L111_[1] + t1*mat.r_L111_[4] + t2*mat.r_L111_[7];
                z = l_new*mat.r_L111_[2] + t1*mat.r_L111_[5] + t2*mat.r_L111_[8];
            }
        }
        else if (v===5) {
            lMax = 2*Math.PI/(LatConst);
            x = state.kx[i] + lMax;
            if (x > lMax) {
                x = -lMax + (x - lMax)
            }
        }
        else if (v===6) {
            lMax = 2*Math.PI/(LatConst);
            y = state.ky[i] + lMax;
            if (y > lMax) {
                y = -lMax + (y - lMax)
            }
        }
    const xPrint = 500/2 + x*kScale;
    const yPrint = 500/2 - y*kScale;
    theContext.beginPath();
    theContext.arc(xPrint, yPrint, 1, 0, 2*Math.PI);
    theContext.fillStyle = '#bdd7e7';
    theContext.fill();
    

    }
    theContext.fillStyle = "white";
    theContext.font= "19px Arial";
    theContext.fillText("k  [100]", 400, 270);
    theContext.font= "13px Arial";
    theContext.fillText("x", 412, 278);   
    theContext.save();                 // save current state
    theContext.fillStyle = "white";
    theContext.font = "19px Arial";
    theContext.translate(270, 100);
    theContext.rotate(-Math.PI / 2);
    theContext.fillText("k  [010]", 0, 0);
    theContext.restore();

    theContext.save();                 // save current state
    theContext.fillStyle = "white";
    theContext.font = "13px Arial";
    theContext.translate(275, 90);
    theContext.rotate(-Math.PI / 2);
    theContext.fillText("y", 0, 0);
    theContext.restore();    
}



function MCLoop(state, scatt) { 
    var eFieldSlider = document.getElementById("eField");
    const Fieldx = Number(eFieldSlider.value)*1e6
    let aveEnHist;
    let totals = {
        EnHist: new Float64Array(NumEnergyLevelStat),
        VelHist: new Float64Array(200),
        En: 0.0,
        Vel: 0.0
    };
    for (let n = 0; n < stepsFrame; n++){
        LoopParticles(state, scatt, -Fieldx)
        totals = statistics(state, totals)
    }
    avgs = averages(totals)
    drawKappa(state, mat)
    drawEnHist(avgs)
    drawVelHist(avgs)
    requestAnimationFrame(() => MCLoop(state, scatt))
}


// // loop
requestAnimationFrame(() => MCLoop(state, scatt));
// // theContext.beginPath();
// // theContext.arc(300, 50, 5s, 0, 2*Math.PI);
// // theContext.fillStyle = "red";
// // theContext.fill();