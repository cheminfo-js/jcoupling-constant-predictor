'use strict';

const OCLE = require('openchemlib-extended');
const getAllCouplings = require('./ocle/getAllCouplings');

class Predictor {
  constructor(options = {db: {}}) {
    let db = options.db;
  }

  predict(molecule, options) {
    // Open the molecule
    // let result = OCLE.Molecule.fromMolfileWithAtomMap(fs.readFileSync(`${structuremol + molid}.mol`).toString());
    // molecule = result.molecule;

    const {
      fromLabel = 'H',
      toLabel = 'H',
      minLength = 1,
      maxLength = 4
    } = options;

    let couplings = getAllCouplings(molecule, { fromLabel, toLabel, minLength, maxLength });

    console.log(couplings)

    /*
    diaIDs = molecule.getGroupedDiastereotopicAtomIDs();
    diaIDs.sort(function (a, b) {
      return b.counter - a.counter;
    });

    for (let k = 0; k < diaIDs.length; k++) {
      diaIDs[k].hose = Util.getHoseCodesFromDiastereotopicID(diaIDs[k].oclID, {
        maxSphereSize: maxSphereSize,
        type: 0
      });
    }

    map = result.map;
    inverseMap = [];
    for (let k = 0; k < map.length; k++) {
      inverseMap.push(map.indexOf(k));
    }*/
  }

  /**
   * Return the euclidean distance between the atom1 and atom2
   * @param {OCL} molecule 
   * @param {Number} atom1 
   * @param {Number} atom2 
   * @returns {Number} 
   */
  distance2(molecule, atom1, atom2) {
    let dx = molecule.getAtomX(atom1) - molecule.getAtomX(atom2);
    let dy = molecule.getAtomY(atom1) - molecule.getAtomY(atom2);
    let dz = molecule.getAtomZ(atom1) - molecule.getAtomZ(atom2);
  
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }
  
  /**
   * Return the dihedral angle between the atoms on the path
   * @param {OCL} molecule 
   * @param {Array} path 
   * @returns {Number} 
   */
  getAngle(molecule, path) {
    let ax = molecule.getAtomX(path[0]) - molecule.getAtomX(path[1]);
    let ay = molecule.getAtomY(path[0]) - molecule.getAtomY(path[1]);
    let az = molecule.getAtomZ(path[0]) - molecule.getAtomZ(path[1]);
  
    let bx = molecule.getAtomX(path[2]) - molecule.getAtomX(path[1]);
    let by = molecule.getAtomY(path[2]) - molecule.getAtomY(path[1]);
    let bz = molecule.getAtomZ(path[2]) - molecule.getAtomZ(path[1]);
  
    return Math.acos((ax * bx + ay * by + az * bz) / ((ax * ax + ay * ay + az * az) * (bx * bx + by * by + bz * bz)));
  }
}

module.exports = Predictor;
