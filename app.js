// DOM Elements
const canvas = document.getElementById('zChart');
const ctx = canvas.getContext('2d');
const mainTitle = document.getElementById('mainTitle');
const distSelect = document.getElementById('distSelect');
const radios = document.getElementsByName('tailTest');
const methodSelect = document.getElementById('methodSelect');
const zMethodInputs = document.getElementById('zMethodInputs');
const pMethodInputs = document.getElementById('pMethodInputs');
const dfSingleContainer = document.getElementById('dfSingleContainer');
const dfFContainer = document.getElementById('dfFContainer');
const feedbackBox = document.getElementById('feedbackBox');
const feedbackText = document.getElementById('feedbackText');
const feedbackExplanation = document.getElementById('feedbackExplanation');
const decisionPrompt = document.getElementById('decisionPrompt');
const lblCrit = document.getElementById('lblCrit');
const lblCalc = document.getElementById('lblCalc');

// Input Elements
const inpDf = document.getElementById('df');
const inpDf1 = document.getElementById('df1');
const inpDf2 = document.getElementById('df2');
const inpStatCrit = document.getElementById('statCrit');
const inpStatCalc = document.getElementById('statCalc');
const inpAlpha = document.getElementById('alpha');
const inpPValue = document.getElementById('pValue');
const saveValuesBtn = document.getElementById('saveValuesBtn');

// Event Listeners
distSelect.addEventListener('change', (e) => {
    const isT = e.target.value === 't';
    const isF = e.target.value === 'f';
    const isC = e.target.value === 'c';

    dfSingleContainer.classList.toggle('hidden', !(isT || isC));
    dfFContainer.classList.toggle('hidden', !isF);

    let statName = 'Z';
    if (isT) statName = 'T';
    if (isF) statName = 'F';
    if (isC) statName = 'χ²';

    mainTitle.innerText = `${statName}-Distribution Hypothesis Testing`;
    lblCrit.innerText = `${statName}-Criteria (${statName.toLowerCase()}_critical)`;
    lblCalc.innerText = `${statName}-Calculation (${statName.toLowerCase()}_calc)`;

    // Update Method dropdown texts
    methodSelect.options[0].innerText = `${statName}-Calculation vs ${statName}-Criteria`;

    // Auto-check right-tail for F and Chi-Square, since it's most common, unless user changed it
    if ((isF || isC) && !Array.from(radios).some(r => r.checked && r.value !== 'two')) {
        document.querySelector('input[name="tailTest"][value="right"]').checked = true;
        updateRadioStyles();
    }

    updatePromptText();
    hideFeedback();
    drawChart();
});

radios.forEach(r => r.addEventListener('change', (e) => {
    updateRadioStyles();
    hideFeedback();
    drawChart();
}));

function updatePromptText() {
    const type = distSelect.value;
    const statName = type === 't' ? 'T' : (type === 'f' ? 'F' : (type === 'c' ? 'χ²' : 'Z'));
    if (methodSelect.value === 'z') {
        decisionPrompt.innerText = `Based on comparing ${statName}-calculation to ${statName}-criteria, what is your conclusion?`;
    } else {
        decisionPrompt.innerText = "Based on comparing P-value to the Level of Significance (α), what is your conclusion?";
    }
}

methodSelect.addEventListener('change', (e) => {
    if (e.target.value === 'z') {
        zMethodInputs.classList.remove('hidden');
        pMethodInputs.classList.add('hidden');
    } else {
        zMethodInputs.classList.add('hidden');
        pMethodInputs.classList.remove('hidden');
    }
    updatePromptText();
    hideFeedback();
    drawChart();
});

saveValuesBtn.addEventListener('click', () => {
    hideFeedback();
    drawChart();
});

function updateRadioStyles() {
    radios.forEach(r => {
        const label = r.parentElement;
        if (r.checked) {
            label.classList.add('bg-blue-50', 'border-blue-200');
        } else {
            label.classList.remove('bg-blue-50', 'border-blue-200');
        }
    });
}

function hideFeedback() {
    feedbackBox.classList.add('hidden');
}

// --- MATH FUNCTIONS ---

// Standard Normal PDF
function normalPDF(x) {
    return (1 / Math.sqrt(2 * Math.PI)) * Math.exp(-0.5 * x * x);
}

// Standard Normal CDF Approximation
function normalCDF(x) {
    let t = 1 / (1 + 0.2316419 * Math.abs(x));
    let d = 0.3989423 * Math.exp(-x * x / 2);
    let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));
    return x > 0 ? 1 - p : p;
}

// Inverse Normal CDF Approximation
function invNorm(p) {
    if (p <= 0.0001) return -4;
    if (p >= 0.9999) return 4;
    const c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
    const d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;
    let t = Math.sqrt(-2 * Math.log(p < 0.5 ? p : 1 - p));
    let z = t - ((c2 * t + c1) * t + c0) / (((d3 * t + d2) * t + d1) * t + 1);
    return p < 0.5 ? -z : z;
}

// Lanczos approximation for log Gamma
function lngamma(z) {
    const cof = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = z;
    let tmp = z + 5.5;
    tmp -= (z + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j <= 5; j++) ser += cof[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / z);
}

function lnbeta(a, b) {
    return lngamma(a) + lngamma(b) - lngamma(a + b);
}

// T-Distribution PDF
function tPDF(t, df) {
    const logNum = lngamma((df + 1) / 2);
    const logDen = lngamma(df / 2) + 0.5 * Math.log(df * Math.PI);
    const logPower = -((df + 1) / 2) * Math.log(1 + (t * t) / df);
    return Math.exp(logNum - logDen + logPower);
}

// T-Distribution CDF
function tCDF(t, df) {
    let isNegative = t < 0;
    t = Math.abs(t);
    let n = 200; // Intervals
    let h = t / n;
    let sum = tPDF(0, df) + tPDF(t, df);
    for (let i = 1; i < n; i++) sum += tPDF(i * h, df) * (i % 2 === 0 ? 2 : 4);
    let areaFromZero = (h / 3) * sum;
    let result = 0.5 + areaFromZero;
    return isNegative ? 1 - result : result;
}

// T-Distribution Inverse CDF
function invT(p, df) {
    if (p <= 0.00001) return -10;
    if (p >= 0.99999) return 10;
    let low = -20, high = 20, mid;
    for (let i = 0; i < 50; i++) {
        mid = (low + high) / 2;
        if (tCDF(mid, df) < p) low = mid;
        else high = mid;
    }
    return mid;
}

// F-Distribution PDF
function fPDF(x, df1, df2) {
    if (x <= 0) return 0;
    const logNum = (df1 / 2) * Math.log(df1) + (df2 / 2) * Math.log(df2) + (df1 / 2 - 1) * Math.log(x);
    const logDen = ((df1 + df2) / 2) * Math.log(df1 * x + df2) + lnbeta(df1 / 2, df2 / 2);
    return Math.exp(logNum - logDen);
}

// F-Distribution CDF (Handles singularities at df1=1 elegantly)
function fCDF(x, df1, df2) {
    if (x <= 0) return 0;
    // Prevent infinity evaluation at x=0 for df1=1 by mapping through T-distribution
    if (df1 === 1) {
        let tVal = Math.sqrt(x);
        return tCDF(tVal, df2) - tCDF(-tVal, df2);
    }
    // Standard numerical integration for df1 >= 2
    let n = 500;
    let start = 1e-6;
    if (x <= start) return 0;
    let h = (x - start) / n;
    let sum = fPDF(start, df1, df2) + fPDF(x, df1, df2);
    for (let i = 1; i < n; i++) sum += fPDF(start + i * h, df1, df2) * (i % 2 === 0 ? 2 : 4);
    return (h / 3) * sum;
}

// F-Distribution Inverse CDF
function invF(p, df1, df2) {
    if (p <= 0.00001) return 0.0001;
    if (p >= 0.99999) return 100;
    let low = 0.0001, high = 100, mid;
    for (let i = 0; i < 60; i++) {
        mid = (low + high) / 2;
        if (fCDF(mid, df1, df2) < p) low = mid;
        else high = mid;
    }
    return mid;
}

// Chi-Square Distribution PDF
function chiPDF(x, df) {
    if (x <= 0) return 0;
    const logNum = (df / 2 - 1) * Math.log(x) - x / 2;
    const logDen = (df / 2) * Math.log(2) + lngamma(df / 2);
    return Math.exp(logNum - logDen);
}

// Chi-Square CDF
function chiCDF(x, df) {
    if (x <= 0) return 0;
    // Prevent infinity evaluation at x=0 for df=1 by mathematically mapping through normal distribution
    if (df === 1) {
        let zVal = Math.sqrt(x);
        return 2 * normalCDF(zVal) - 1;
    }
    // Standard numerical integration for df >= 2
    let n = 500;
    let start = 1e-6;
    if (x <= start) return 0;
    let h = (x - start) / n;
    let sum = chiPDF(start, df) + chiPDF(x, df);
    for (let i = 1; i < n; i++) sum += chiPDF(start + i * h, df) * (i % 2 === 0 ? 2 : 4);
    return (h / 3) * sum;
}

// Chi-Square Inverse CDF
function invChi(p, df) {
    if (p <= 0.00001) return 0.0001;
    if (p >= 0.99999) return df * 5 + 20; // safe dynamic upper limit
    let low = 0.0001, high = Math.max(100, df * 5), mid;
    for (let i = 0; i < 60; i++) {
        mid = (low + high) / 2;
        if (chiCDF(mid, df) < p) low = mid;
        else high = mid;
    }
    return mid;
}

// --- DRAWING FUNCTION ---

function drawChart() {
    const w = canvas.width;
    const h = canvas.height;
    const padding = 40;

    // Get Current Values
    const distType = distSelect.value;
    const testType = Array.from(radios).find(r => r.checked).value;
    const method = methodSelect.value;
    const df = parseFloat(inpDf.value) || 10;
    const df1 = parseFloat(inpDf1.value) || 5;
    const df2 = parseFloat(inpDf2.value) || 10;

    // Choose appropriate math functions
    let pdfFn, invFn;
    if (distType === 'f') {
        pdfFn = (x) => fPDF(x, df1, df2);
        invFn = (p) => invF(p, df1, df2);
    } else if (distType === 'c') {
        pdfFn = (x) => chiPDF(x, df);
        invFn = (p) => invChi(p, df);
    } else if (distType === 't') {
        pdfFn = (x) => tPDF(x, df);
        invFn = (p) => invT(p, df);
    } else {
        pdfFn = normalPDF;
        invFn = invNorm;
    }

    // Define X axis ranges based on distribution
    let xMin = -4;
    let xMax = 4;
    if (distType === 'f') {
        xMin = 0;
        xMax = 6;
    } else if (distType === 'c') {
        xMin = 0;
        xMax = Math.ceil(df + 4 * Math.sqrt(2 * df)); // Dynamically scale up to Mean + 4 SD
        if (xMax < 10) xMax = 10;
    }
    const xRange = xMax - xMin;

    // Coordinate mapping limits
    const mapX = (val) => padding + ((val - xMin) / xRange) * (w - 2 * padding);

    let statCrit = 0;
    let statCalc = 0;

    if (method === 'z') {
        statCrit = parseFloat(inpStatCrit.value) || 0;
        if (distType !== 'f' && distType !== 'c') statCrit = Math.abs(statCrit); // Z/T are symmetric inputs usually
        statCalc = parseFloat(inpStatCalc.value) || 0;
    } else {
        const alpha = parseFloat(inpAlpha.value) || 0.05;
        const pVal = parseFloat(inpPValue.value) || 0.5;

        // Map P-values to corresponding scores for accurate visual plotting
        if (testType === 'left') {
            statCrit = (distType === 'f' || distType === 'c') ? invFn(alpha) : Math.abs(invFn(alpha));
            statCalc = invFn(pVal);
        } else if (testType === 'right') {
            statCrit = (distType === 'f' || distType === 'c') ? invFn(1 - alpha) : Math.abs(invFn(1 - alpha));
            statCalc = invFn(1 - pVal);
        } else {
            statCrit = (distType === 'f' || distType === 'c') ? invFn(1 - alpha / 2) : Math.abs(invFn(1 - alpha / 2));
            statCalc = invFn(1 - pVal / 2); // Plotting positive side representation
        }
    }

    // Determine maximum Y height to fit the curve nicely
    let maxPdf = 0;
    const curvePoints = [];
    const step = xRange / 300; // 300 responsive points for smoothness across dynamic ranges
    for (let i = xMin; i <= xMax; i += step) {
        let evalX = i;
        if ((distType === 'f' || distType === 'c') && evalX === 0) evalX = 0.001; // Avoid strict 0
        let y = pdfFn(evalX);

        // Cap visual height for extreme peaks around 0 for df=1 (F and Chi-Sq)
        if (distType === 'f' && df1 === 1 && y > 1.5) y = 1.5;
        if (distType === 'c' && df === 1 && y > 1.5) y = 1.5;

        if (y > maxPdf && y < Infinity) maxPdf = y;
        curvePoints.push({ x: evalX, y: y });
    }
    maxPdf = Math.max(maxPdf * 1.1, 0.1); // Add headroom
    if (distType === 'f' && df1 === 1) maxPdf = 1.0;
    if (distType === 'c' && df === 1) maxPdf = 1.0;

    const mapY = (val) => h - padding - (Math.min(val, maxPdf) / maxPdf) * (h - 2 * padding);

    ctx.clearRect(0, 0, w, h);

    // Draw Base X-Axis
    // Draw Standard Ticks (Dynamically spaced based on xRange to prevent overlap)
    ctx.fillStyle = '#64748b';
    ctx.font = '12px sans-serif';
    ctx.textAlign = 'center';

    let tickStep = 1;
    if (xRange > 15) tickStep = 5;
    if (xRange > 40) tickStep = 10;
    if (xRange > 80) tickStep = 20;

    for (let i = Math.ceil(xMin / tickStep) * tickStep; i <= Math.floor(xMax); i += tickStep) {
        const xPos = mapX(i);
        ctx.beginPath();
        ctx.moveTo(xPos, h - padding);
        ctx.lineTo(xPos, h - padding + 5);
        ctx.stroke();
        ctx.fillText(i, xPos, h - padding + 20);
    }

    // Helper to shade regions seamlessly
    function shadeRegion(start, end) {
        if (start > xMax || end < xMin) return;
        start = Math.max(start, xMin);
        end = Math.min(end, xMax);

        ctx.beginPath();
        ctx.moveTo(mapX(start), mapY(0));
        for (let x = start; x <= end; x += step) {
            let evalX = (distType === 'f' || distType === 'c') && x === 0 ? 0.001 : x;
            let y = pdfFn(evalX);
            if (distType === 'f' && df1 === 1 && y > 1.5) y = 1.5;
            if (distType === 'c' && df === 1 && y > 1.5) y = 1.5;
            ctx.lineTo(mapX(x), mapY(y));
        }
        let endY = pdfFn(end);
        if (distType === 'f' && df1 === 1 && endY > 1.5) endY = 1.5;
        if (distType === 'c' && df === 1 && endY > 1.5) endY = 1.5;
        ctx.lineTo(mapX(end), mapY(endY));
        ctx.lineTo(mapX(end), mapY(0));
        ctx.closePath();
        ctx.fillStyle = 'rgba(239, 68, 68, 0.4)'; // Tailwind red-500 with opacity
        ctx.fill();
    }

    // Determine regions to shade
    if (distType === 'f' || distType === 'c') {
        if (testType === 'left') {
            shadeRegion(xMin, statCrit);
        } else if (testType === 'right') {
            shadeRegion(statCrit, xMax);
        } else if (testType === 'two') {
            let lowerBound = statCrit;
            if (method === 'z') {
                // Infer the missing lower/upper tail based on user's single critical input
                let pValCrit = distType === 'f' ? fCDF(statCrit, df1, df2) : chiCDF(statCrit, df);
                if (pValCrit >= 0.5) {
                    lowerBound = distType === 'f' ? invF(1 - pValCrit, df1, df2) : invChi(1 - pValCrit, df);
                    shadeRegion(xMin, lowerBound);
                    shadeRegion(statCrit, xMax);
                } else {
                    let upperBound = distType === 'f' ? invF(1 - pValCrit, df1, df2) : invChi(1 - pValCrit, df);
                    shadeRegion(xMin, statCrit);
                    shadeRegion(upperBound, xMax);
                }
            } else {
                // P-method accurately provides both mappings automatically
                shadeRegion(xMin, invFn((parseFloat(inpAlpha.value) || 0.05) / 2));
                shadeRegion(invFn(1 - (parseFloat(inpAlpha.value) || 0.05) / 2), xMax);
            }
        }
    } else {
        // Z and T Distributions
        if (testType === 'left') {
            shadeRegion(xMin, -statCrit);
        } else if (testType === 'right') {
            shadeRegion(statCrit, xMax);
        } else if (testType === 'two') {
            shadeRegion(xMin, -statCrit);
            shadeRegion(statCrit, xMax);
        }
    }

    // Draw the Curve Outline
    ctx.beginPath();
    ctx.moveTo(mapX(curvePoints[0].x), mapY(curvePoints[0].y));
    curvePoints.forEach(p => {
        ctx.lineTo(mapX(p.x), mapY(p.y));
    });
    ctx.strokeStyle = '#1e293b'; // slate-800
    ctx.lineWidth = 2.5;
    ctx.stroke();

    // Draw Test Statistic
    const zx = mapX(statCalc);
    if (zx >= mapX(xMin) && zx <= mapX(xMax)) {
        // Draw dashed line
        ctx.beginPath();
        ctx.moveTo(zx, h - padding);
        ctx.lineTo(zx, mapY(Math.min(pdfFn(statCalc), maxPdf)));
        ctx.strokeStyle = '#2563eb'; // blue-600
        ctx.lineWidth = 2;
        ctx.setLineDash([5, 5]);
        ctx.stroke();
        ctx.setLineDash([]);

        // Draw dot
        ctx.beginPath();
        ctx.arc(zx, mapY(Math.min(pdfFn(statCalc), maxPdf)), 4, 0, 2 * Math.PI);
        ctx.fillStyle = '#2563eb';
        ctx.fill();

        // Label
        ctx.fillStyle = '#1d4ed8';
        ctx.font = 'bold 12px sans-serif';
        // Adjust text position so it doesn't clip
        const alignLeft = zx > w / 2;
        ctx.textAlign = alignLeft ? 'right' : 'left';
        const xOffset = alignLeft ? -8 : 8;
        const statPrefix = distType === 't' ? 'T' : (distType === 'f' ? 'F' : (distType === 'c' ? 'χ²' : 'Z'));
        const labelText = method === 'z' ? `${statPrefix}-calc (${statCalc.toFixed(2)})` : `P-value equivalent`;
        ctx.fillText(labelText, zx + xOffset, mapY(Math.min(pdfFn(statCalc), maxPdf)) - 10);
    }
}

// --- DECISION CHECKER ---

function checkDecision(userChoice) {
    const distType = distSelect.value;
    const statName = distType === 't' ? 'T' : (distType === 'f' ? 'F' : (distType === 'c' ? 'χ²' : 'Z'));
    const testType = Array.from(radios).find(r => r.checked).value;
    const method = methodSelect.value;

    let isReject = false;
    let explanation = '';

    if (method === 'z') {
        const statCrit = parseFloat(inpStatCrit.value) || 0;
        const statCalc = parseFloat(inpStatCalc.value) || 0;

        if (distType === 'f' || distType === 'c') {
            const df1 = parseFloat(inpDf1.value) || 5;
            const df2 = parseFloat(inpDf2.value) || 10;
            const dfC = parseFloat(inpDf.value) || 10;

            if (testType === 'left') {
                isReject = statCalc <= statCrit;
                explanation = `For a left-tailed ${statName}-test, we reject if ${statName}-calc ( ${statCalc} ) ≤ ${statName}-crit ( ${statCrit} ).`;
            } else if (testType === 'right') {
                isReject = statCalc >= statCrit;
                explanation = `For a right-tailed ${statName}-test, we reject if ${statName}-calc ( ${statCalc} ) ≥ ${statName}-crit ( ${statCrit} ).`;
            } else {
                // Two-tailed asymmetric distribution logic
                let pValCrit = distType === 'f' ? fCDF(statCrit, df1, df2) : chiCDF(statCrit, dfC);
                let upper, lower;
                if (pValCrit >= 0.5) {
                    upper = statCrit;
                    lower = distType === 'f' ? invF(1 - pValCrit, df1, df2) : invChi(1 - pValCrit, dfC);
                } else {
                    lower = statCrit;
                    upper = distType === 'f' ? invF(1 - pValCrit, df1, df2) : invChi(1 - pValCrit, dfC);
                }
                isReject = statCalc >= upper || statCalc <= lower;
                explanation = `For a two-tailed ${statName}-test, we check both extremes. Reject if ${statName}-calc ( ${statCalc} ) ≥ ${statName}-upper ( ${upper.toFixed(3)} ) OR ${statName}-calc ≤ ${statName}-lower ( ${lower.toFixed(3)} ).`;
            }
        } else {
            const absCrit = Math.abs(statCrit);
            if (testType === 'left') {
                isReject = statCalc <= -absCrit;
                explanation = `For a left-tailed test, we reject if ${statName}-calc ( ${statCalc} ) ≤ -${statName}-crit ( -${absCrit} ).`;
            } else if (testType === 'right') {
                isReject = statCalc >= absCrit;
                explanation = `For a right-tailed test, we reject if ${statName}-calc ( ${statCalc} ) ≥ ${statName}-crit ( ${absCrit} ).`;
            } else {
                isReject = Math.abs(statCalc) >= absCrit;
                explanation = `For a two-tailed test, we reject if |${statName}-calc| ( ${Math.abs(statCalc)} ) ≥ ${statName}-crit ( ${absCrit} ).`;
            }
        }

    } else {
        const alpha = parseFloat(inpAlpha.value) || 0;
        const pVal = parseFloat(inpPValue.value) || 0;

        // P-value method is universal across all distributions!
        isReject = pVal <= alpha;
        explanation = `Using the P-value method, we reject H₀ if P-value ( ${pVal} ) ≤ Level of Significance α ( ${alpha} ).`;
    }

    const correctChoice = isReject ? 'reject' : 'fail';
    const isCorrect = userChoice === correctChoice;

    // Display Feedback UI
    feedbackBox.classList.remove('hidden', 'bg-green-100', 'text-green-800', 'bg-red-100', 'text-red-800');

    if (isCorrect) {
        feedbackBox.classList.add('bg-green-100', 'text-green-800');
        feedbackText.innerHTML = '✅ <strong>Correct Decision!</strong>';
    } else {
        feedbackBox.classList.add('bg-red-100', 'text-red-800');
        feedbackText.innerHTML = '❌ <strong>Incorrect Decision.</strong>';
    }

    feedbackExplanation.innerText = explanation + (isReject ? " Therefore, we Reject the Null Hypothesis." : " Therefore, we Fail to Reject the Null Hypothesis.");
}

// Initialize on load
window.onload = () => {
    updateRadioStyles();
    drawChart();
};