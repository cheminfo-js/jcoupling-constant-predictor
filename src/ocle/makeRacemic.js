'use strict';

const { Molecule } = require('openchemlib-extended');

function makeRacemic(molecule) {
  // if we don't calculate this we have 2 epimers
  molecule.ensureHelperArrays(Molecule.cHelperParities);
  // we need to make one group "AND" for chiral (to force to racemic, this means diastereotopic and not enantiotopic)
  for (let i = 0; i < molecule.getAllAtoms(); i++) {
    if (molecule.getAtomParity(i) !== Molecule.cAtomParityNone) {
      molecule.setAtomESR(i, Molecule.cESRTypeAnd, 0); // changed to group 0; TLS 9.Nov.2015
    }
  }
}

module.exports = makeRacemic;