'use strict'

const OCLE = require('openchemlib-extended');
const nmrPredictor = require('nmr-predictor');
const Predictor = require('./index');
const fs = require('fs');
const path = require('path');
const stat = require('ml-stat/array');
const histogram = require('./histogram');

let db = JSON.parse(fs.readFileSync('data/cheminfo-absHH.json').toString());

let folder = '/Users/acastillo//Documents/dataNMR/spinus/'
let data = fs.readdirSync(folder).filter(file => file.indexOf('.mol') >= 0);

let nMols = data.length;
let nSamples = 200;
let p = new Predictor({ db });

// let stats = { min: Number.MAX_VALUE, max: Number.MIN_VALUE, sum: 0, count: 0 };
let stats = [];


for (let n = 0; n < nSamples; n++) {
 // console.log(n)
    let randomSample = Math.round(Math.random() * nMols);
    let molfile = fs.readFileSync(path.join(folder, data[randomSample])).toString();
    // console.log(molfile)
    let spinus = JSON.parse(fs.readFileSync(path.join(folder, data[randomSample].replace('.mol', '.json'))));
    let molmap = OCLE.Molecule.fromMolfileWithAtomMap(molfile);
    // console.log(JSON.stringify(spinus))

    let hydrogens = molmap.molecule.getGroupedDiastereotopicAtomIDs({ atomLabel: 'H' });
    hydrogens = hydrogens.reduce((acum, value) => {
        acum.push(...value.atoms);
        return acum;
    }, []);

    // console.log(hydrogens);
    let jSpinus = symmetrizeInPlace(asMatrix(spinus, hydrogens));
    
    let couplings = p.predict3D(molmap.molecule, { mapper: x => x }).filter(x => x.fromDiaID !== x.toDiaID);
    let jCheminfo = symmetrizeInPlace(p.asMatrix(couplings, hydrogens, x => Math.abs(x.median)  * 1.5));


    //console.log(JSON.stringify(jCheminfo))
    //console.log(JSON.stringify(jSpinus))

    let result = compare(jCheminfo, jSpinus);
    stats.push(...result);
    //stats.sum += result.sum;
    //stats.count += result.count;
    //console.log(result);
    //console.log(result.sum / result.count);
    //console.log(result.absError + ' ' + result.count);
}

let resume = getStats(stats);
fs.writeFileSync("stats.txt", stats.join());

let hist =  histogram({data: stats, bins: linspace(resume.min, resume.max, 100)});
console.log(resume.median + ' ' +resume.mean);
console.log(hist.map(x => x.y).join(','));



function compare(a, b) {
    let lng = a.length;
    /*let sum = 0;
    let countA = 0;
    let countB = 0;
    let countAB = 0;*/
    let diff = [];
    for (let i = 0; i < lng; i++) {
        for (let j = i + 1; j < lng; j++) {
            if(a[i][j]) {
                if (b[i][j]) {
                    diff.push(Math.abs(a[i][j] - b[i][j]));
                    //sum += Math.abs(a[i][j] - b[i][j]);
                    //countAB++;
                } else {
                    diff.push(a[i][j]);
                    //sum += Math.abs(a[i][j]);
                    //countA++;
                }
            } else if (b[i][j]) {
                diff.push(b[i][j]);
                //sum += Math.abs(b[i][j]);
                //countB++;
            }
        }
    }
    return diff;
    //return {sum, countAB, countA, countB, count: countAB + countA + countB };

}

function asMatrix(couplings, atomIDs) {
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
                if (ji.assignment[0] * 1 === atom2) {
                    return  Math.abs(ji.coupling);
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
            if (couplings[from][to] !== undefined && couplings[to][from] !==undefined) {
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
  