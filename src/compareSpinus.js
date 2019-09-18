'use strict'

const OCLE = require('openchemlib-extended');
const nmrPredictor = require('nmr-predictor');
const Predictor = require('./index');
const fs = require('fs');
const path = require('path');


let db = JSON.parse(fs.readFileSync('data/cheminfoHH.json').toString());

let folder = '/Users/acastillo//Documents/dataNMR/spinus/'
let data = fs.readdirSync(folder).filter(file => file.indexOf('.mol') >= 0);

let nMols = data.length;
let nSamples = 1;
let p = new Predictor({ db });

let stats = { min: Number.MAX_VALUE, max: Number.MIN_VALUE, sum: 0 };

for (let n = 0; n < nSamples; n++) {
    let randomSample = n;//Math.round(Math.random() * nMols);
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
    let jCheminfo = symmetrizeInPlace(p.asMatrix(couplings, hydrogens, x => Math.abs(x.mean)));


    console.log(JSON.stringify(jCheminfo))
    // console.log(jSpinus)

}

console.log(stats);

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
