const OCLE = require('openchemlib-extended');
const fs = require('fs');

const Predictor = require('../src/index');

let result = OCLE.Molecule.fromMolfileWithAtomMap(fs.readFileSync('data/ethylbenzene3d.mol').toString());
let db = JSON.parse(fs.readFileSync('data/cheminfoHH.json').toString());

let p = new Predictor({ db });
let couplings = p.predict3D(result.molecule, { mapper: x => x });
let hydrogens = result.molecule.getGroupedDiastereotopicAtomIDs({ atomLabel: 'H' });
hydrogens = hydrogens.reduce((acum, value) => {
    acum.push(...value.atoms);
    return acum;
}, []);
let jcc = p.asMatrix(couplings, hydrogens, x => Math.abs(x.mean));

// console.log(JSON.stringify(couplings));
console.log(jcc);