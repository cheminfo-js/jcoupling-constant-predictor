const nmrPredictor = require('nmr-predictor'),
const jccPredictor = require('./index');
const fs = require('fs');
const path = require('path');

let db = JSON.parse(fs.readFileSync('data/cheminfoHH.json').toString());

let folder = '~/Documents/dataNMR/spinus/'
let data = fs.readdirSync(folder).filter(file => file.indexOf('.mol')>=0);

let nMols = data.length;

let nSamples = 100;
let stats = {min: Number.MAX_VALUE, max: Number.MIN_VALUE, sum: 0}
for(let n = 0; n < nSamples; n++) {
    let randomSample = Math.round(Math.random() * nMols);
    let molfile = fs.readFileSync(path.join(folder, data[randomSample])).toString();
    let spinus = JSON.parse(fs.readFileSync(path.join(folder, data[randomSample].replace('.mol', '.json'))));

    let molmap = OCLE.Molecule.fromMolfileWithAtomMap(molfile);

}

console.log(stats);
