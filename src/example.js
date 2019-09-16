const OCLE = require('openchemlib-extended');
const fs = require('fs');

const Predictor = require('../src/index');


let result = OCLE.Molecule.fromMolfileWithAtomMap(fs.readFileSync('data/ethylbenzene3d.mol').toString());

let p = new Predictor();

p.predict(result.molecule, {});