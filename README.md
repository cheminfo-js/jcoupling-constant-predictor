# j-coupling-constant-predictor

[![NPM version][npm-image]][npm-url]
[![build status][travis-image]][travis-url]
[![Test coverage][codecov-image]][codecov-url]
[![David deps][david-image]][david-url]
[![npm download][download-image]][download-url]

Predictor of 1,2,3,4 and 5 J H-H coupling constants

## Installation

`$ npm i j-coupling-constant-predictor`

## [API Documentation](https://cheminfo-js.github.io/j-coupling-constant-predictor/)

## Example

```js
const predictor = require('j-coupling-constant-predictor');
const fs = require('fs');
const OCLE = require('openchemlib-extended');

let molfile = fs.readFileSync('moleculeWithExpandedHydrogens.mol').toString();
let molmap = OCLE.Molecule.fromMolfileWithAtomMap(molfile);

//Predict the coupling constants using the 3D information of the molecule. Use the mean of the most similar
//entries as value for the coupling constant. You can use median aswell.
let couplings = predictor.predict3D(molmap.molecule, {maxLength: 6, mapper: x => x.mean });

// Remove the couplings between chemically equivalent atoms
couplings = couplings.filter(x => x.fromDiaID !== x.toDiaID);

```

## License

[MIT](./LICENSE)

[npm-image]: https://img.shields.io/npm/v/j-coupling-constant-predictor.svg?style=flat-square
[npm-url]: https://www.npmjs.com/package/j-coupling-constant-predictor
[travis-image]: https://img.shields.io/travis/com/cheminfo-js/j-coupling-constant-predictor/master.svg?style=flat-square
[travis-url]: https://travis-ci.com/cheminfo-js/j-coupling-constant-predictor
[codecov-image]: https://img.shields.io/codecov/c/github/cheminfo-js/j-coupling-constant-predictor.svg?style=flat-square
[codecov-url]: https://codecov.io/gh/cheminfo-js/j-coupling-constant-predictor
[david-image]: https://img.shields.io/david/cheminfo-js/j-coupling-constant-predictor.svg?style=flat-square
[david-url]: https://david-dm.org/cheminfo-js/j-coupling-constant-predictor
[download-image]: https://img.shields.io/npm/dm/j-coupling-constant-predictor.svg?style=flat-square
[download-url]: https://www.npmjs.com/package/j-coupling-constant-predictor
