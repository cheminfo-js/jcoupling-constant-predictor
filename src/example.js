const OCLE = require('openchemlib-extended');
const fs = require('fs');

const Predictor = require('../src/index');


let result = OCLE.Molecule.fromMolfileWithAtomMap(fs.readFileSync('data/ethylbenzene3d.mol').toString());
let db = JSON.parse(fs.readFileSync('data/cheminfoHH.json').toString());
//console.log(db);

let p = new Predictor({db});

p.predict3D(result.molecule, {});