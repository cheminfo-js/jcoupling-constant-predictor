'use strict';

const OCLE = require('openchemlib-extended');
const getAllCouplings = require('./ocle/getAllCouplings');
const Util = OCLE.Util;

class Predictor {
  constructor(options = { db: {} }) {
    this.db = options.db;
    this.maxSphereSize = 5;
  }

  predictPlain(molecule, options) {

  }

  predict3D(molecule, options) {
    const {
      fromLabel = 'H',
      toLabel = 'H',
      minLength = 2,
      maxLength = 3
    } = options;

    let couplings = getAllCouplings(molecule, { fromLabel, toLabel, minLength, maxLength });
    // console.log(couplings);


    let diaIDs = molecule.getGroupedDiastereotopicAtomIDs();
    diaIDs.sort(function (a, b) {
      return b.counter - a.counter;
    });

    // HOSE codes from Diasterotopic atom IDs.
    for (let k = 0; k < diaIDs.length; k++) {
      diaIDs[k].hose = Util.getHoseCodesFromDiastereotopicID(diaIDs[k].oclID, {
        maxSphereSize: this.maxSphereSize,
        type: 0
      });
    }

    // console.log(diaIDs);
    /* map = result.map;
      inverseMap = [];
      for (let k = 0; k < map.length; k++) {
        inverseMap.push(map.indexOf(k));
      }*/

    couplings.forEach(chemPair => {
      // molecule.getPath(atoms, inverseMap[example[2]], inverseMap[example[3]], 3);
      let type = `${chemPair.pathLength}J${chemPair.fromLabel}${chemPair.toLabel}`;
      let dbt = this.db[type];
      if (dbt) {
        let diaID = diaIDs.find(x => x.oclID == chemPair.fromDiaID);
        let hose = diaID.hose;

        // Huggly query
        let pred = null;
        if (hose[this.maxSphereSize - 1]) {
          pred = dbt[2][hose[this.maxSphereSize - 1]];
          if (!pred) {
            pred = dbt[1][hose[this.maxSphereSize - 2]];
            if (!pred) {
              pred = dbt[0][hose[this.maxSphereSize - 3]];
              if (pred)
                pred.lvl = this.maxSphereSize - 2; 
            } else {
              pred.lvl = this.maxSphereSize - 1; 
            }
          } else {
            pred.lvl = this.maxSphereSize; 
          }
        }
        pred.couplings = [];
        chemPair.fromTo.forEach(magPair => {
          let atoms = [];
          molecule.getPath(atoms, magPair[0], magPair[1], 3);
          let feature = 0;
          console.log(atoms);
          if (atoms.length === 4) {
            console.log('Torsion');
            feature = molecule.calculateTorsion(atoms);
            pred.couplings.push(pred.cop2[0] + pred.cop2[1] * Math.cos(feature) + pred.cop2[2] * Math.pow(Math.cos(feature), 2));
          } else {
            console.log('Distance');
            feature = this.distance2(molecule, atoms[0], atoms[2]);
            pred.couplings.push(pred.cop2[0] + pred.cop2[1] * feature + pred.cop2[2] * Math.pow(feature, 2));
          }
        });
        console.log(type);
        console.log(pred);
      }
    });

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
   * Return the angle between the atoms on the path
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
