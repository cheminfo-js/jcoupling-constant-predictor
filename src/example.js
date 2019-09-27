const OCLE = require('openchemlib-extended');
const fs = require('fs');

const Predictor = require('../src/index');

let result = OCLE.Molecule.fromMolfileWithAtomMap(fs.readFileSync('data/ethylbenzene3d.mol').toString());
let db = JSON.parse(fs.readFileSync('data/cheminfo-abs-spinusHH.json').toString());

let p = new Predictor({ db });
let couplings = p.predict3D(result.molecule, { mapper: x => x }).filter(x => (x.fromDiaID !== x.toDiaID) && x.j);;
let diaIDsH = result.molecule.getGroupedDiastereotopicAtomIDs({ atomLabel: 'H' });
console.log(diaIDsH);
let hydrogens = diaIDsH.reduce((acum, value) => {
    acum.push(...value.atoms);
    return acum;
}, []);
console.log(hydrogens);
let jCheminfo = p.asMatrix(couplings, hydrogens, crazyFit);


// console.log(JSON.stringify(couplings));
console.log(jCheminfo);

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