'use strict';

const OCLE = require('openchemlib-extended');

const xAtomicNumber = OCLE.Molecule.getAtomicNoFromLabel('X');

function changeAtom(molecule, iAtom) {
  molecule.setAtomCustomLabel(iAtom, `${molecule.getAtomLabel(iAtom)}*`);
  if (molecule.getAtomicNo(iAtom) === 1) {
    molecule.setAtomicNo(iAtom, xAtomicNumber);
  } else {
    // we can not use X because we would have problems with valencies if it is
    // expanded hydrogens or not
    // we can not only use a custom label because it does not count for the canonisation
    molecule.setAtomMass(iAtom, molecule.getAtomMass(iAtom) + 5);
  }
}

module.exports = changeAtom;