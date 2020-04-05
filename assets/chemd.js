var canvas = document.getElementById("chemdoodle-viewer");
let transformUnitCell = new ChemDoodle.TransformCanvas3D('transformUnitCell', 200, 200);
// set up styles
transformUnitCell.styles.set3DRepresentation('Ball and Stick');
transformUnitCell.styles.backgroundColor = 'black';
transformUnitCell.styles.atoms_displayLabels_3D = true;
transformUnitCell.styles.shapes_color = '#fff';
// create a gold atom
let mol = new ChemDoodle.structures.Molecule();
mol.atoms.push(new ChemDoodle.structures.Atom('Au'));
// create a cubic unit cell
let unitCell = new ChemDoodle.structures.d3.UnitCell({ o: [-.5, -.5, -.5], x: [.5, -.5, -.5], y: [-.5, .5, -.5], z: [-.5, -.5, .5], xy: [.5, .5, -.5], xz: [.5, -.5, .5], yz: [-.5, .5, .5], xyz: [.5, .5, .5] });
// add the objects to the scene
transformUnitCell.loadContent([mol], [unitCell]);