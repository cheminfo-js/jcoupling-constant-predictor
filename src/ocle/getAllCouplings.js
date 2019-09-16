'use strict';

var floydWarshall = require('ml-floyd-warshall');
var Matrix = require('ml-matrix').Matrix;

const OCLE = require('openchemlib-extended');
const changeAtom = require('./changeAtom');
const makeRacemic = require('./makeRacemic');
const getAllPaths = require('./getAllPaths');

/**
 * Returns an array of all the different atom diaIDs that are connected
 * {OCLE} molecule
 * {object} [options={}]
 * {string} [options.fromLabel='H']
 * {string} [options.toLabel='H']
 * {number} [options.minLength=1]
 * {number} [options.maxLength=4]
 */

function getAllCouplings(molecule, options = {}) {
  const {
    fromLabel = 'H',
    toLabel = 'H',
    minLength = 1,
    maxLength = 4
  } = options;
  let paths =getAllPaths(molecule, {
    fromLabel,
    toLabel,
    minLength,
    maxLength
  });

  const minSphereSize = 0;
  const maxSphereSize = 2;

  const connectivityMatrix = molecule.getConnectivityMatrix();
  const pathLengthMatrix = floydWarshall(new Matrix(connectivityMatrix));

  let fragment = new OCLE.Molecule(0, 0);
  for (let path of paths) { // .slice(0, 1)) {
    path.info = [];
    for (let fromTo of path.fromTo) {
      let atoms = [];
      molecule.getPath(atoms, fromTo[0], fromTo[1], path.pathLength);
      let torsion;
      if (atoms.length === 4) {
        torsion = molecule.calculateTorsion(atoms);
      }
      path.info.push({
        atoms,
        torsion
      });

      if (!path.code) {
        path.code = [];
        let tmpMolecule = molecule.getCompactCopy();
        makeRacemic(tmpMolecule);
        changeAtom(tmpMolecule, atoms[0]);
        changeAtom(tmpMolecule, atoms[atoms.length - 1]);
        let atomMask = new Array(tmpMolecule.getAllAtoms()).fill(false);
        for (let atom of atoms) {
          atomMask[atom] = true;
        }
        let previousSphere = '';
        for (var sphere = 0; sphere <= maxSphereSize; sphere++) {
          for (let atom of atoms) {
            let distances = pathLengthMatrix[atom];
            for (let i = 0; i < distances.length; i++) {
              if (distances[i] === sphere && !atomMask[i]) {
                if (tmpMolecule.getAtomLabel(i) !== 'H' ||
                 (tmpMolecule.getAtomLabel(i) === 'H' && tmpMolecule.getAtomLabel(tmpMolecule.getConnAtom(i, 0)) !== 'C')) {
                  atomMask[i] = true;
                }
              }
            }
          }
          if (previousSphere !== atomMask.join()) {
            previousSphere = atomMask.join();
            tmpMolecule.copyMoleculeByAtoms(fragment, atomMask, true, null);
            if (sphere >= minSphereSize) {
              path.code.push(
                fragment.getCanonizedIDCode(
                  OCLE.Molecule.CANONIZER_ENCODE_ATOM_CUSTOM_LABELS
                )
              );
            }
          }
        }
      }
    }
  }
  return paths;
}

module.exports = getAllCouplings;
