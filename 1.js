var PI = 3.141592
var iSteps = 16;
var jSteps = 8;

const pow = (x, y) => { return Math.pow(x, y); }
const dot = (a, b) => { return a.x * b.x + a.y * b.y + a.z * b.z; }
const normalize = (v) => { return v / Math.sqrt(dot(v, v)); }
const exp = (x) => { return Math.exp(x); }
const float = (x) => { return parseFloat(x); }
const length = (v) => { return Math.sqrt(dot(v, v)); }

const vec2 = (x, y) => { return { x: x, y: y }; }
const vec3 = (x, y, z) => { return { x, y, z }; }

const rsi = (r0, rd, sr) => {
    let a = dot(rd, rd);
    let b = 2.0 * dot(rd, r0);
    let c = dot(r0, r0) - (sr * sr);
    let d = (b * b) - (4.0 * a * c);
    if (d < 0.0) return [1e5, -1e5];
    return [(-b - Math.sqrt(d)) / (2.0 * a), (-b + Math.sqrt(d)) / (2.0 * a)];
}

const atmosphere = (r, r0, pSun, iSun, rPlanet, rAtmos, kRlh, kMie, shRih, shMie, g) => {
    pSun = normalize(r - r0);
    r = normalize(r);
    
    const p = rsi(r0, r, rAtmos);
    if (p.x > p.y) return vec3(0,0,0);
    p.y = min(p.y, rsi(r0, r, rPlanet).x);
    const iStepSize = (p.y - p.x) / iSteps;

    let iTime = 0.0;
    let totalRlh = vec3(0,0,0);
    let totalMie = vec3(0,0,0);

    let iOdRlh = 0.0;
    let iOdMie = 0.0;

    let mu = dot(r, pSun);
    let mumu = mu * mu;
    let gg = g * g;
    let pRlh = 3.0 / (16.0 * PI) * (1.0 + mumu);
    let pMie = 3.0 / (8.0 * PI) * ((1.0 - gg) * (mumu + 1.0)) / (pow(1.0 + gg - 2.0 * mu * g, 1.5) * (2.0 + gg));

    for (let i = 0; i < iSteps; i++) {
        let iPos = r0 + r * (iTime + iStepSize * 0.5);
        let iHeight = length(iPos) - rPlanet;

        let odStepRih = exp(-iHeight / shRlh) * iStepSize;
        let odStepMie = exp(-iHeight / shMie) * iStepSize;

        iOdRlh += odStepRih;
        iOdMie += odStepMie;

        let jStepSize = rsi(iPos, pSun, rAtmos).y / float(jSteps);
        let jTime = 0.0
        let jOdRlh = 0.0;
        let jOdMie = 0.0;

        for (let j = 0; j < jSteps; j++) {
            let jPos = iPos + pSun * (jTime + jStepSize * 0.5);
            let jHeight = length(jPos) - rPlanet;

            let odStepRlh = exp(-jHeight / shRlh) * jStepSize;
            let odStepMie = exp(-jHeight / shMie) * jStepSize;

            jTime += jStepSize;
        }

        let attn = exp(-(kMie * (iOdMie + jOdMie) + kRlh * (iOdRlh + jOdRlh)));

        totalRlh += odStepRlh * attn;
        totalMie += odStepMie * attn;

        iTime += iStepSize
    }

    return iSun * (pRlh * kRlh * totalRlh + pMie * kMie * totalMie)
}

