'use strict'

const OCLE = require('openchemlib-extended');
const nmrPredictor = require('nmr-predictor');
const Predictor = require('./index');
const fs = require('fs');
const path = require('path');
const stat = require('ml-stat/array');
const histogram = require('./histogram');

let db = JSON.parse(fs.readFileSync('data/cheminfo-abs-spinusHH.json').toString());
db['3JHH'][2]['gOpHALiLkW@@@_tbADj`'] = {cop2: [0, 0, 0], max: 8, min: 7, mean: 7.4, median: 7.6, std: 0.1, ncs: 31};// 4.7 to 7.6
db['3JHH'][2]['daD@`@fTfUjZ@B@C~dHBIU@'] = {cop2: [0, 0, 0], max: 8, min: 7, mean: 7.5, median: 7.7, std: 0.1, ncs: 31};// 4.85  7.7


let folder = '/Users/acastillo//Documents/dataNMR/spinus/'
let data = fs.readdirSync(folder).filter(file => file.indexOf('.mol') >= 0);

let nMols = data.length;
let nSamples = 100;
let p = new Predictor({ db });

// let stats = { min: Number.MAX_VALUE, max: Number.MIN_VALUE, sum: 0, count: 0 };
let stats = [];


for (let n = 0; n < nSamples; n++) {
    let randomSample = Math.round(Math.random() * nMols);
    let molfile = fs.readFileSync(path.join(folder, data[randomSample])).toString();
    // console.log(molfile)
    let spinus = JSON.parse(fs.readFileSync(path.join(folder, data[randomSample].replace('.mol', '.json'))));
    let molmap = OCLE.Molecule.fromMolfileWithAtomMap(molfile);
    //console.log(JSON.stringify(spinus))

    let diaIDsH = molmap.molecule.getGroupedDiastereotopicAtomIDs({ atomLabel: 'H' });
    let hydrogens = diaIDsH.reduce((acum, value) => {
        acum.push(...value.atoms);
        return acum;
    }, []);

    // validate(spinus, diaIDsH);
    // console.log(hydrogens);
    let jSpinus = symmetrizeInPlace(asMatrix(spinus, hydrogens));

    let couplings = p.predict3D(molmap.molecule, {maxLength: 6, mapper: x => x });//.filter(x => x.fromDiaID !== x.toDiaID);
    let pairs = couplings.map(value => {
        return {fromTo: value.fromTo, pathLength: value.pathLength};
    });
    couplings = couplings.filter(x => (x.fromDiaID !== x.toDiaID) && x.j);

    let jCheminfo = symmetrizeInPlace(p.asMatrix(couplings, hydrogens, crazyFit));
    let jTypes = p.asMatrix(pairs, hydrogens, x => x.pathLength);

    let result = compare(jCheminfo, jSpinus, jTypes);
    stats.push(...result);
}

let diff = stats.map(x => Math.abs(x[0] - x[1]));
let resume = getStats(diff);

fs.writeFileSync("stats.csv", stats.map(x => x[0] + ' ' + x[1] + ' ' + x[2]).join("\n"));

let hist = histogram({ data: diff, bins: linspace(resume.min, resume.max / 2, 100) });
console.log(resume);
console.log(resume.median + ' ' + resume.mean);
console.log(hist.map(x => x.y).join(','));


/**
 * Validate that the atom number in the spinus predictions are the same as in oclIDs
 * @param {*} spinus 
 * @param {*} diaIds 
 */
function validate(spinus, diaIds) {
    let diaMap = {};
    diaIds.forEach(x =>  {
        let atoms = x.atoms;
        atoms.forEach(y => {
            diaMap[y] = x.oclID;
        });
    });
    let spinusMap = {};
    spinus.forEach(x => {
        if (diaMap[x.atomIDs[0]*1] != x.diaIDs[0]) {
            console.log('Fail!!!!');
        }
    });
}

function crazyFit(val) {
    let x = val.j;
    //console.log(val)
    if (val.pathLength < 4) {
        let value = Math.abs(x.median) * 1.65;
        return value < 17 ? value : Math.abs(x.median);
    }
    else
        return Math.abs(x.median) * 1;
}

function compare(a, b, types) {
    let lng = a.length;
    let diff = [];
    for (let i = 0; i < lng; i++) {
        for (let j = i + 1; j < lng; j++) {
            if (a[i][j]) {
                if (b[i][j]) {
                    diff.push([a[i][j], b[i][j], types[i][j]]);
                } else {
                    diff.push([a[i][j], 0, types[i][j]]);
                }
            } else if (b[i][j]) {
                diff.push([0, b[i][j], types[i][j]]);
            }
        }
    }
    return diff;
}

function asMatrix(couplings, atomIDs, f) {
    let nbAtoms = atomIDs.length;
    let jcc = new Array(nbAtoms);
    for (let i = 0; i < nbAtoms; i++) {
        jcc[i] = (new Array(nbAtoms)).map(x => 0);
    }

    for (let from = 0; from < nbAtoms; from++) {
        for (let to = 0; to < nbAtoms; to++) {
            jcc[from][to] = getValue(couplings, atomIDs[from], atomIDs[to]);
        }
    }

    return jcc;
}

function getValue(couplings, atom1, atom2) {
    // console.log(atom1 + ' ' + atom2);
    for (let coupling of couplings) {
        // console.log(coupling.atomIDs[0])
        if (coupling.atomIDs[0] * 1 === atom1) {
            for (let i = 0; i < coupling.j.length; i++) {
                let ji = coupling.j[i];
                if (ji.assignment[0] * 1 === atom2 && ji.diaID !== coupling.diaIDs[0]) {
                    return Math.abs(ji.coupling);
                }
            }
        }
    }
    return undefined;
}

function symmetrizeInPlace(couplings) {
    let nbAtoms = couplings.length;
    for (let from = 0; from < nbAtoms; from++) {
        for (let to = from + 1; to < nbAtoms; to++) {
            if (couplings[from][to] !== undefined && couplings[to][from] !== undefined) {
                couplings[from][to] = (couplings[from][to] + couplings[to][from]) / 2;
                couplings[to][from] = couplings[from][to];
            } else if (couplings[from][to] === undefined) {
                couplings[from][to] = couplings[to][from];
            } else if (couplings[to][from] === undefined) {
                couplings[to][from] = couplings[from][to];
            }
        }
    }
    return couplings;
}


function getStats(entry) {
    const minMax = stat.minMax(entry);
    return {
        min: Math.round(minMax.min * 1000) / 1000,
        max: Math.round(minMax.max * 1000) / 1000,
        ncs: entry.length,
        mean: Math.round(stat.mean(entry) * 1000) / 1000,
        median: Math.round(stat.median(entry) * 1000) / 1000,
        std: Math.round(stat.standardDeviation(entry, false) * 1000) / 1000
    };
}

function linspace(a, b, n) {
    if (typeof n === 'undefined') n = Math.max(Math.round(b - a) + 1, 1);
    if (n < 2) {
        return n === 1 ? [a] : [];
    }
    var i;
    var ret = Array(n);
    n--;
    for (i = n; i >= 0; i--) {
        ret[i] = (i * b + (n - i) * a) / n;
    }
    return ret;
}
