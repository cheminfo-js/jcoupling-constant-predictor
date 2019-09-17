'use strict';

const OCLE = require('openchemlib-extended');
const getAllCouplings = require('./ocle/getAllCouplings');
const Util = OCLE.Util;

class Predictor {
  constructor(options = { db: {} }) {
    this.db = options.db;
    this.maxSphereSize = 5;
  }

  /**
   * Predict the JHH coupling constants for the input molecule. It uses a database
   * compiled from the data of ab-inition calculations provided for the kaggle competition
   * https://www.kaggle.com/c/champs-scalar-coupling and a set of rules for 4JHH extracted from 
   * https://www.chem.wisc.edu/areas/reich/nmr/05-hmr-06-4j.htm
   * The functions uses the 3D geometry if it is available, but it is the user who must decide
   * if it can be used or not. It could differ from the 3D prediction, But I don't see how right now-
   * @param {OCLE} molecule 
   * @param {object} options 
   */
  predictPlain(molecule, options) {
    return this.predict3D(molecule, options);
  }

  /**
   * Predict the JHH coupling constants for the input molecule. It uses a database
   * compiled from the data of ab-inition calculations provided for the kaggle competition
   * https://www.kaggle.com/c/champs-scalar-coupling and a set of rules for 4JHH extracted from 
   * https://www.chem.wisc.edu/areas/reich/nmr/05-hmr-06-4j.htm
   * The functions uses the 3D geometry if it is available, but it is the user who must decide
   * if it can be used or not. 
   * @param {OCLE} molecule 
   * @param {object} options 
   */
  predict3D(molecule, options) {
    const {
      fromLabel = 'H',
      toLabel = 'H',
      minLength = 2,
      maxLength = 4
    } = options;

    let couplings = getAllCouplings(molecule, { fromLabel, toLabel, minLength, maxLength });

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
      if (type == '4JHH') {
        let pred = {couplings: []};
        chemPair.j = pred;
        chemPair.fromTo.forEach(magPair => {
          let atoms = [];
          molecule.getPath(atoms, magPair[0], magPair[1], 4);
          if (!this.isAttachedToHeteroAtom(molecule, magPair[0]) && !this.isAttachedToHeteroAtom(molecule, magPair[1])) {
            this.predict4JHH(molecule, atoms, pred);
          }
        });
      } else {
        if (this.db[type]) {
          let pred = this.query(this.db[type], diaIDs.find(x => x.oclID == chemPair.fromDiaID).hose, this.maxSphereSize);
          pred.couplings = [];
          chemPair.j = pred;
          chemPair.fromTo.forEach(magPair => {
            let atoms = [];
            molecule.getPath(atoms, magPair[0], magPair[1], 3);
            let feature = 0;
            switch (atoms.length) {
              case 4:
                feature = Math.cos(molecule.calculateTorsion(atoms));
                break;
              case 3:
                feature = this.distance2(molecule, atoms[0], atoms[2]);
                break;
              case 2:
                feature = this.distance2(molecule, atoms[0], atoms[1]);
                break;
            }
            pred.couplings.push(pred.cop2[0] + pred.cop2[1] * feature + pred.cop2[2] * Math.pow(feature, 2));
          });
        }
      }
    });
    return couplings;
  }

  /**
   * Give a value for 4JHH based on: https://www.chem.wisc.edu/areas/reich/nmr/05-hmr-06-4j.htm
   * @param {OCLE} molecule 
   * @param {Array} atoms 
   * @param {object} pred 
   */
  predict4JHH(molecule, atoms, pred) {
    // TODO: W couplings
    // Allylic-coupling
    if (this.isAllylic(molecule, atoms)) {
      // Betweeen -3 +3 Hz. Return the mean of the absolut value 1.5
      pred.couplings.push(1.5);
      pred.kind = "allylic";
    }
    // Propargylic
    if (this.isPropargylic(molecule, atoms)) {
      // Between +2 +4. Return the mean of the absolut value 3
      pred.couplings.push(3);
      pred.kind = "propargylic";
    }
    // Allenic
    if (this.isAllenic(molecule, atoms)) {
      // Between +6 +7. Return the mean of the absolut value 6.5
      pred.couplings.push(6.5);
      pred.kind = "allenic";
    }
    // Meta-coupling in aromatic compounds
    if (this.isMetaAromatic(molecule, atoms)) {
      // Between +1 +3. Return the mean of the absolut value 2
      pred.couplings.push(2);
      pred.kind = "mata-aromatic";
    }
  }

  /**
   * Query the JCC db with the given hose code
   * @param {Array} dbt 
   * @param {Array} hose 
   * @param {Number} maxSphereSize 
   */
  query(dbt, hose, maxSphereSize) {
    let pred = null;
    if (hose[maxSphereSize - 1]) {
      pred = dbt[2][hose[maxSphereSize - 1]];
      if (!pred) {
        pred = dbt[1][hose[maxSphereSize - 2]];
        if (!pred) {
          pred = dbt[0][hose[maxSphereSize - 3]];
          if (pred)
            pred.lvl = maxSphereSize - 2;
        } else {
          pred.lvl = maxSphereSize - 1;
        }
      } else {
        pred.lvl = maxSphereSize;
      }
    }
    return pred;
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
   * Return the angle between the atoms on the 3 length-path
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

  isSingleBond(molecule, atom1, atom2) {
    var bond = molecule.getBond(atom1, atom2);
    var bondType = molecule.getBondType(bond);
    return bondType === 1;
  }

  isDoubleBond(molecule, atom1, atom2) {
    var bond = molecule.getBond(atom1, atom2);
    var bondType = molecule.getBondType(bond);
    return bondType === 2;
  }

  isTripleBond(molecule, atom1, atom2) {
    var bond = molecule.getBond(atom1, atom2);
    var bondType = molecule.getBondType(bond);
    return bondType === 4;
  }

  isDoubleOrTripleBond(molecule, atom1, atom2) {
    var bond = molecule.getBond(atom1, atom2);
    var bondType = molecule.getBondType(bond);
    return bondType === 2 || bondType === 4;
  }

  isAromatic(molecule, atom1, atom2) {
    var bond = molecule.getBond(atom1, atom2);
    return molecule.isAromaticBond(bond);
  }

  isMetaAromatic(molecule, atoms) {
    if (this.isDoubleBond(molecule, atoms[1], atoms[2]) &&
      this.isSingleBond(molecule, atoms[2], atoms[3]) &&
      this.isAromatic(molecule, atoms[1], atoms[2])) {
      return true;
    }
    if (this.isDoubleBond(molecule, atoms[2], atoms[3]) &&
      this.isSingleBond(molecule, atoms[1], atoms[2]) &&
      this.isAromatic(molecule, atoms[2], atoms[3])) {
      return true;
    }
    return false;
  }

  isAllylic(molecule, atoms) {
    if (this.isDoubleBond(molecule, atoms[1], atoms[2]) &&
      this.isSingleBond(molecule, atoms[2], atoms[3]) &&
      !this.isAromatic(molecule, atoms[1], atoms[2])) {
      return true;
    }
    if (this.isDoubleBond(molecule, atoms[2], atoms[3]) &&
      this.isSingleBond(molecule, atoms[1], atoms[2]) &&
      !this.isAromatic(molecule, atoms[2], atoms[3])) {
      return true;
    }
    return false;
  }

  isPropargylic(molecule, atoms) {
    if (this.isTripleBond(molecule, atoms[1], atoms[2]) &&
      this.isSingleBond(molecule, atoms[2], atoms[3]) &&
      !this.isAromatic(molecule, atoms[1], atoms[2])) {
      return true;
    }
    if (this.isTripleBond(molecule, atoms[2], atoms[3]) &&
      this.isSingleBond(molecule, atoms[1], atoms[2]) &&
      !this.isAromatic(molecule, atoms[2], atoms[3])) {
      return true;
    }
    return false;
  }

  isAllenic(molecule, atoms) {
    if (this.isDoubleBond(molecule, atoms[1], atoms[2]) &&
      this.isDoubleBond(molecule, atoms[2], atoms[3]) &&
      !this.isAromatic(molecule, atoms[1], atoms[2])) {
      return true;
    }
    if (this.isDoubleBond(molecule, atoms[2], atoms[3]) &&
      this.isDoubleBond(molecule, atoms[1], atoms[2]) &&
      !this.isAromatic(molecule, atoms[2], atoms[3])) {
      return true;
    }
    return false;
  }

  isAttachedToHeteroAtom(molecule, atom) {
    var nbConnectedAtoms = molecule.getAllConnAtoms(atom);
    for (var j = 0; j < nbConnectedAtoms; j++) {
      var connAtom = molecule.getConnAtom(atom, j);
      if (!(molecule.getAtomLabel(connAtom) === 'C')) {
        return true;
      }
    }
    return false;
  }

}

module.exports = Predictor;
