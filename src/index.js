'use strict';

const OCLE = require('openchemlib-extended');
const getAllCouplings = require('./ocle/getAllCouplings');
const Util = OCLE.Util;

class Predictor {
  /**
   * Constructor
   * @param {object} options Actually it is not an option. You must provide a db
   */
  constructor(options = { db: {} }) {
    this.db = options.db;
    // HOSE 5
    // gOpHALiLkW@@@_tbADj`    4.7 to 7.6
    // daD@`@fTfUjZ@B@C~dHBIU@ 4.85  7.7
    //daF@`FBYRYVkh@`@OzP`HeT  5
    this.maxSphereSize = 5;
  }

  /**
   * Convert the array of couplings objects into a numeric matrix
   * containing all the coupling constant values between each pair of 
   * atomIDs. The matrix could be asymmetric
   * @param {Array} couplings 
   * @param {Array} atomIDs 
   * @param {function} mapper 
   * @returns {array} A numeric matrix with coupling constants
   */
  asMatrix(couplings, atomIDs, mapper) {
    let nbAtoms = atomIDs.length;
    if (!mapper)
      mapper = j => j;
    let jcc = new Array(nbAtoms);
    for (let i = 0; i < nbAtoms; i++) {
      jcc[i] = (new Array(nbAtoms)).map(x => 0);
    }

    for (let from = 0; from < nbAtoms; from++) {
      for (let to = 0; to < nbAtoms; to++) {
        jcc[from][to] = getValue(couplings, atomIDs[from], atomIDs[to], mapper);
      }
    }

    return jcc;
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
   * compiled from the data of ab-initio calculations provided for the kaggle competition
   * https://www.kaggle.com/c/champs-scalar-coupling and a set of rules for 4JHH extracted from 
   * https://www.chem.wisc.edu/areas/reich/nmr/05-hmr-06-4j.htm
   * The functions uses the 3D geometry if it is available, but it is the user who must decide
   * if it can be used or not. 
   * @param {OCLE} molecule 
   * @param {object} options 
   * @options {string} fromLabel
   * @options {string} toLabel
   * @options {number} minLength
   * @options {number} maxLength 
   * @options {function} mapper A function to transform out the j component of the prediction.  
   */
  predict3D(molecule, options) {
    const {
      fromLabel = 'H',
      toLabel = 'H',
      minLength = 2,
      maxLength = 4,
      mapper = x => x
    } = options;

    let couplings = getAllCouplings(molecule, { fromLabel, toLabel, minLength, maxLength });

    let diaIDs = molecule.getGroupedDiastereotopicAtomIDs();
    diaIDs.sort(function (a, b) {
      return b.counter - a.counter;
    });

    // HOSE codes from Diasterotopic atom IDs.
    let diaIDsMap = {};
    for (let k = 0; k < diaIDs.length; k++) {
      diaIDs[k].hose = Util.getHoseCodesFromDiastereotopicID(diaIDs[k].oclID, {
        maxSphereSize: this.maxSphereSize,
        type: 0
      });
      diaIDsMap[diaIDs[k].oclID] = diaIDs[k];
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
      let pred = {};
      /*if (type == '4JHH' || type == '5JHH') {
        chemPair.fromTo.forEach(magPair => {
          let atoms = [];
          molecule.getPath(atoms, magPair[0], magPair[1], 5);
          if (!isAttachedToHeteroAtom(molecule, magPair[0]) && !isAttachedToHeteroAtom(molecule, magPair[1])) {
            if (type == '4JHH')
              this.predict4JHH(molecule, atoms, pred);
            if (type == '5JHH')
              this.predict5JHH(molecule, atoms, pred);
          }
        });
      } else */{
        if (this.db[type]) {
          pred = this.query(this.db, type, diaIDsMap[chemPair.fromDiaID].hose, diaIDsMap[chemPair.toDiaID].hose, this.maxSphereSize);
          if (pred) {
            pred.reg = [];
            chemPair.fromTo.forEach(magPair => {
              let atoms = [];
              molecule.getPath(atoms, magPair[0], magPair[1], 3);
              let feature = 0;
              switch (atoms.length) {
                case 4:
                  feature = Math.cos(molecule.calculateTorsion(atoms));
                  break;
                case 3:
                  feature = distance2(molecule, atoms[0], atoms[2]);
                  break;
                case 2:
                  feature = distance2(molecule, atoms[0], atoms[1]);
                  break;
              }
              pred.reg.push(pred.cop2[0] + pred.cop2[1] * feature + pred.cop2[2] * Math.pow(feature, 2));
            });
          } else {
            if (!pred || !pred.lvl) {
              if (type === '4JHH' || type === '5JHH') {
                pred = {};
                chemPair.fromTo.forEach(magPair => {
                  let atoms = [];
                  molecule.getPath(atoms, magPair[0], magPair[1], 4);
                  if (!isAttachedToHeteroAtom(molecule, magPair[0]) && !isAttachedToHeteroAtom(molecule, magPair[1])) {
                    if (type == '4JHH' )
                      pred = this.predict4JHH(molecule, atoms, pred);
                    if (type == '5JHH' )
                      pred = this.predict5JHH(molecule, atoms, pred);
                  }
                });
              }
            }
          }
        }
      }
      if (pred && Object.keys(pred).length !== 0)
        chemPair.j = mapper(pred);
    });
    return couplings;// couplings.filter(entry => entry.j);
  }

  /**
   * Give a value for 4JHH based on: https://www.chem.wisc.edu/areas/reich/nmr/05-hmr-06-4j.htm
   * @param {OCLE} molecule 
   * @param {Array} atoms 
   * @param {object} pred 
   */
  predict4JHH(molecule, atoms, pred) {
    // Meta-coupling in aromatic compounds
    // Disabled because it was producing a lot missmatches with spinus
    if (isMetaAromatic(molecule, atoms)) {
      // Between +1 +3. Return the mean of the absolut value 2
      pred.mean = 2;
      pred.median = 1.6;
      pred.min = 1;
      pred.max = 3;
      pred.kind = "mata-aromatic";
    }
    // TODO: W couplings
    // Allylic-coupling
    /*if (isAllylic(molecule, atoms)) {
      // Betweeen -3 +3 Hz. Return the mean of the absolut value 1.5
      pred.mean = 1.5;
      pred.median = 1.5;
      pred.min = -3;
      pred.max = 3;
      pred.kind = "allylic";
    }
    // Propargylic
    if (isPropargylic(molecule, atoms)) {
      // Between +2 +4. Return the mean of the absolut value 3
      pred.mean = 3;
      pred.median = 3;
      pred.min = 2;
      pred.max = 4;
      pred.kind = "propargylic";
    }
    // Allenic
    if (isAllenic(molecule, atoms)) {
      // Between +6 +7. Return the mean of the absolut value 6.5
      pred.mean = 6.5;
      pred.median = 6.5;
      pred.min = 6;
      pred.max = 7;
      pred.kind = "allenic";
    }*/
    return pred;
  }

  /**
 * Give a value for 5JHH based on: https://www.chem.wisc.edu/areas/reich/nmr/05-hmr-06-4j.htm
 * @param {OCLE} molecule 
 * @param {Array} atoms 
 * @param {object} pred 
 */
  predict5JHH(molecule, atoms, pred) {
    // console.log(atoms)
    if (isHomoAllylic(molecule, atoms)) {
      // Between 0 +8. Return 0.5 because I think
      pred.mean = 4;
      pred.median = 0,5;
      pred.min = 0;
      pred.max = 8;
      pred.lvl = 0;
      pred.kind = "homoallylic";
    }
    if (isHomoPropargylic(molecule, atoms)) {
      // Between +2 +4. Return the mean of the absolut value 3
      pred.mean = 3;
      pred.median = 1.5;
      pred.min = 2;
      pred.max = 4;
      pred.kind = "homopropargylic";
    }
  }

  /**
   * Query the JCC db with the given hose code
   * @param {Array} dbt 
   * @param {Array} hose 
   * @param {Number} maxSphereSize 
   */
  query(db, type, hoseFrom, hoseTo, maxSphereSize) {
    let dbt = db[type];
    let pred = null;
    let lng = type.substring(0, 1) * 1;
    if (lng < 4) {
      if (hoseFrom[maxSphereSize - 1]) {
        pred = dbt[2][hoseFrom[maxSphereSize - 1]];
        if (!pred) {
          pred = dbt[1][hoseFrom[maxSphereSize - 2]];
          if (!pred) {
            pred = dbt[0][hoseFrom[maxSphereSize - 3]];
            if (pred)
              pred = Object.assign({lvl: maxSphereSize - 2}, pred);
          } else {
            pred = Object.assign({lvl: maxSphereSize - 1}, pred);
          }
        } else {
          pred = Object.assign({lvl: maxSphereSize}, pred);
        }
      }
    } else {
      if (hoseFrom[maxSphereSize - 1] && hoseTo[maxSphereSize - 1]) {
        let key = canCat(hoseFrom[maxSphereSize - 1], hoseTo[maxSphereSize - 1]);
        pred = dbt[2][key];
        if (pred)
          pred = Object.assign({lvl: maxSphereSize}, pred);
      }
    }
    /*if (pred && pred.lvl)
      return pred;
    else {
      return {mean: 2.9, median: 2.925, min: 2.925, max: 4.43, lvl: 0, cop2: [2.925, 0, 0]}
    } */
    return pred;
  }
}


/**
 * Cat the strings a and b doing min(a, b) + max(a, b)
 * @param {*} a 
 * @param {*} b 
 */
function canCat(a, b) {
  if (a < b)
    return a + b;
  else
    return b + a;
}

/**
 * Get the value of the coupling constant betweeen atom1 and atom2. The mapper
 * specify how to convert the object to the output value
 * @param {Array} couplings 
 * @param {number} atom1 
 * @param {number} atom2 
 * @param {function} mapper 
 * @returns {*} 
 */
function getValue(couplings, atom1, atom2, mapper) {
  for (let coupling of couplings) {
    for (let pair of coupling.fromTo) {
      if (pair[0] == atom1 && pair[1] == atom2) {
        return mapper(coupling);
      }
    }
  }
}

/**
 * Return the euclidean distance between the atom1 and atom2
 * @param {OCL} molecule 
 * @param {Number} atom1 
 * @param {Number} atom2 
 * @returns {Number} 
 */
function distance2(molecule, atom1, atom2) {
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
function getAngle(molecule, path) {
  let ax = molecule.getAtomX(path[0]) - molecule.getAtomX(path[1]);
  let ay = molecule.getAtomY(path[0]) - molecule.getAtomY(path[1]);
  let az = molecule.getAtomZ(path[0]) - molecule.getAtomZ(path[1]);

  let bx = molecule.getAtomX(path[2]) - molecule.getAtomX(path[1]);
  let by = molecule.getAtomY(path[2]) - molecule.getAtomY(path[1]);
  let bz = molecule.getAtomZ(path[2]) - molecule.getAtomZ(path[1]);

  return Math.acos((ax * bx + ay * by + az * bz) / ((ax * ax + ay * ay + az * az) * (bx * bx + by * by + bz * bz)));
}

function isSingleBond(molecule, atom1, atom2) {
  var bond = molecule.getBond(atom1, atom2);
  var bondType = molecule.getBondType(bond);
  return bondType === 1;
}

function isDoubleBond(molecule, atom1, atom2) {
  var bond = molecule.getBond(atom1, atom2);
  var bondType = molecule.getBondType(bond);
  return bondType === 2;
}

function isTripleBond(molecule, atom1, atom2) {
  var bond = molecule.getBond(atom1, atom2);
  var bondType = molecule.getBondType(bond);
  return bondType === 4;
}

function isDoubleOrTripleBond(molecule, atom1, atom2) {
  var bond = molecule.getBond(atom1, atom2);
  var bondType = molecule.getBondType(bond);
  return bondType === 2 || bondType === 4;
}

function isAromatic(molecule, atom1, atom2) {
  var bond = molecule.getBond(atom1, atom2);
  return molecule.isAromaticBond(bond);
}

function isMetaAromatic(molecule, atoms) {
  if (isAromatic(molecule, atoms[1], atoms[2]) &&
      isAromatic(molecule, atoms[2], atoms[3])) {
        if (isDoubleBond(molecule, atoms[1], atoms[2]) &&
          isSingleBond(molecule, atoms[2], atoms[3])) {
          return true;
        }
        if (isDoubleBond(molecule, atoms[2], atoms[3]) &&
          isSingleBond(molecule, atoms[1], atoms[2])) {
          return true;
        }
  }
  return false;
}

function isAllylic(molecule, atoms) {
  if (isDoubleBond(molecule, atoms[1], atoms[2]) &&
    isSingleBond(molecule, atoms[2], atoms[3]) &&
    !isAromatic(molecule, atoms[1], atoms[2])) {
    return true;
  }
  if (isDoubleBond(molecule, atoms[2], atoms[3]) &&
    isSingleBond(molecule, atoms[1], atoms[2]) &&
    !isAromatic(molecule, atoms[2], atoms[3])) {
    return true;
  }
  return false;
}

function isPropargylic(molecule, atoms) {
  if (isTripleBond(molecule, atoms[1], atoms[2]) &&
    isSingleBond(molecule, atoms[2], atoms[3]) &&
    !isAromatic(molecule, atoms[1], atoms[2])) {
    return true;
  }
  if (isTripleBond(molecule, atoms[2], atoms[3]) &&
    isSingleBond(molecule, atoms[1], atoms[2]) &&
    !isAromatic(molecule, atoms[2], atoms[3])) {
    return true;
  }
  return false;
}

function isAllenic(molecule, atoms) {
  if (isDoubleBond(molecule, atoms[1], atoms[2]) &&
    isDoubleBond(molecule, atoms[2], atoms[3]) &&
    !isAromatic(molecule, atoms[1], atoms[2])) {
    return true;
  }
  if (isDoubleBond(molecule, atoms[2], atoms[3]) &&
    isDoubleBond(molecule, atoms[1], atoms[2]) &&
    !isAromatic(molecule, atoms[2], atoms[3])) {
    return true;
  }
  return false;
}

function isAttachedToHeteroAtom(molecule, atom) {
  var nbConnectedAtoms = molecule.getAllConnAtoms(atom);
  for (var j = 0; j < nbConnectedAtoms; j++) {
    var connAtom = molecule.getConnAtom(atom, j);
    if (!(molecule.getAtomLabel(connAtom) === 'C')) {
      return true;
    }
  }
  return false;
}

function isHomoAllylic(molecule, atoms) {
  if (isSingleBond(molecule, atoms[1], atoms[2]) &&
    isSingleBond(molecule, atoms[3], atoms[4]) &&
    isDoubleBond(molecule, atoms[2], atoms[3])) {
    return true;
  }
  return false;
}

function isHomoPropargylic(molecule, atoms) {
  if (isSingleBond(molecule, atoms[1], atoms[2]) &&
    isSingleBond(molecule, atoms[3], atoms[4]) &&
    isTripleBond(molecule, atoms[2], atoms[3])) {
    return true;
  }
  return false;
}

module.exports = Predictor;

