'use strict';

var floydWarshall = require('ml-floyd-warshall');
var Matrix = require('ml-matrix').Matrix;

/**
 * Returns an array of all the different atom diaIDs that are connected
 *
 * This function exports a list of all shortest paths among each pair of atoms in the molecule.
 * It is possible to filter out the list by atomLabel and min/max length. The resulting array will
 * contains both, the diasterotopicAtomID and the atomIDs between the atom in the path. This is important
 * because it allows to differentiate chemically equivalent atoms from magentically equivalent atoms.
 * {object} [options={}]
 * {string} [options.fromLabel='']:q!
 * {string} [options.toLabel='']
 * {number} [options.minLength=1]
 * {number} [options.maxLength=4]
*/

module.exports = function getAllPaths(molecule, options = {}) {
  const {
    fromLabel = '',
    toLabel = '',
    minLength = 1,
    maxLength = 4
  } = options;

  // we need to find all the atoms 'fromLabel' and 'toLabel'
  var results = {};
  var diaIDs = molecule.getDiastereotopicAtomIDs();

  var connectivityMatrix = molecule.getConnectivityMatrix();
  // TODO have a package that allows to convert the connectivityMatrix to a distanceMatrix
  var pathLengthMatrix = floydWarshall(new Matrix(connectivityMatrix));

  for (var from = 0; from < molecule.getAllAtoms(); from++) {
    for (var to = 0; to < molecule.getAllAtoms(); to++) {
      if (!fromLabel || molecule.getAtomLabel(from) === fromLabel) {
        if (!toLabel || molecule.getAtomLabel(to) === toLabel) {
          var pathLength = pathLengthMatrix[from][to];
          var key = `${diaIDs[from]}_${diaIDs[to]}_${pathLength}`;
          if (pathLength >= minLength && pathLength <= maxLength) {
            if (!results[key]) {
              results[key] = {
                fromDiaID: diaIDs[from],
                toDiaID: diaIDs[to],
                fromTo: [[from, to]],
                fromLabel: molecule.getAtomLabel(from),
                toLabel: molecule.getAtomLabel(to),
                pathLength: pathLength
              };
            } else {
              results[key].fromTo.push([from, to]);
            }
          }
        }
      }
    }
  }

  var finalResults = [];
  for (let key in results) {
    finalResults.push(results[key]);
  }
  return finalResults;
};