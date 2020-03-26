import os
import ast
import argparse

from opencmiss.zinc.context import Context
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from opencmiss.zinc.status import OK as ZINC_OK


class ProgramArguments(object):
    pass


def elmult(u, v):
    return [(u[i] * v[i]) for i in range(len(u))]


def scale_coordinates(field, scale, nodeset='nodes', isotropic=True):
    number_of_components = field.getNumberOfComponents()
    if (number_of_components != 2) and (number_of_components != 3):
        print('field has invalid number of components')
        return False
    if len(scale) != number_of_components:
        print('invalid matrix number of columns or offset size')
        return False
    if field.getCoordinateSystemType() != Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN:
        print('field is not rectangular cartesian')
        return False
    fe_field = field.castFiniteElement()
    if not fe_field.isValid():
        print('field is not finite element field type')
        return False
    success = True
    fm = field.getFieldmodule()
    fm.beginChange()
    cache = fm.createFieldcache()
    if nodeset == 'nodes':
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    else:
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
    node_template = nodes.createNodetemplate()
    node_iter = nodes.createNodeiterator()
    node = node_iter.next()
    if isotropic:
        sum_scale = scale[0] + scale[1] + scale[2]
        mean_scale = sum_scale / 3
        final_scale = [mean_scale] * 3
    else:
        final_scale = scale
    while node.isValid():
        node_template.defineFieldFromNode(fe_field, node)
        cache.setNode(node)
        for derivative in [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                           Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3,
                           Node.VALUE_LABEL_D2_DS2DS3, Node.VALUE_LABEL_D3_DS1DS2DS3]:
            versions = node_template.getValueNumberOfVersions(fe_field, -1, derivative)
            for v in range(1, versions + 1):
                result, values = fe_field.getNodeParameters(cache, -1, derivative, v, number_of_components)
                if result != ZINC_OK:
                    success = False
                else:
                    new_values = elmult(final_scale, values)
                    result = fe_field.setNodeParameters(cache, -1, derivative, v, new_values)
                    if result != ZINC_OK:
                        success = False
        node = node_iter.next()
    fm.endChange()
    if not success:
        print('zincutils.transformCoordinates: failed to get/set some values')
    return success


class TransformExGeometry(object):

    def __init__(self, filename, output_filename, scale, marker_type='nodes'):
        self._context = Context("TrasformMesh")
        self._region = self._context.getDefaultRegion()
        self._region.readFile(filename)

        self._output_file = output_filename

        self._marker_type = marker_type

        self._scale = scale
        self._scale_geometry()

    def _get_mesh(self):
        fm = self._region.getFieldmodule()
        for dimension in range(3, 0, -1):
            mesh = fm.findMeshByDimension(dimension)
            if mesh.getSize() > 0:
                return mesh
        raise ValueError('Model contains no mesh')

    def _scale_geometry(self):
        fm = self._region.getFieldmodule()
        scale_coordinates(fm.findFieldByName("coordinates"),
                          self._scale,
                          nodeset='nodes',
                          isotropic=True)

        scale_coordinates(fm.findFieldByName("marker_coordinates"),
                          self._scale,
                          nodeset=self._marker_type,
                          isotropic=True)
        self._write()

    def _get_model_coordinate_field(self):
        mesh = self._get_mesh()
        element = mesh.createElementiterator().next()
        if not element.isValid():
            raise ValueError('Model contains no elements')
        fm = self._region.getFieldmodule()
        cache = fm.createFieldcache()
        cache.setElement(element)
        field_iter = fm.createFielditerator()
        field = field_iter.next()
        while field.isValid():
            if field.isTypeCoordinate() and (field.getNumberOfComponents() <= 3):
                if field.isDefinedAtLocation(cache):
                    return field
            field = field_iter.next()
        raise ValueError('Could not determine model coordinate field')

    def _write(self):
        self._region.writeFile(self._output_file)


def main():
    args = parse_args()
    if os.path.exists(args.input_ex):
        if args.output_ex is None:
            filename = os.path.basename(args.input_ex)
            dirname = os.path.dirname(args.input_ex)
            output_ex = os.path.join(dirname, filename.split('.')[0] + '_scaled.' + filename.split('.')[1])
        else:
            output_ex = args.output_ex

    if args.scale is None:
        return

    scale = ast.literal_eval(args.scale)

    if args.marker_type is None:
        mtype = 'nodes'
    else:
        mtype = args.marker_type

    teg = TransformExGeometry(args.input_ex, output_ex, scale, mtype)


def parse_args():
    parser = argparse.ArgumentParser(description="Transform (currently only scale) EX data.")
    parser.add_argument("input_ex", help="Location of the input EX file.")
    parser.add_argument("-o", "--output-ex", help="Location of the output ex file. "
                                            "[defaults to the location of the input file if not set.]")
    parser.add_argument("-s", "--scale", help="A list to scale the EX e.g. [0.1, 0.1, 0.1]"
                                        "Default is None.")
    parser.add_argument("-m", "--marker_type", help="The domain of the marker coordinates e.g. nodes or datapoint."
                                        "Default is nodes.")

    program_arguments = ProgramArguments()
    parser.parse_args(namespace=program_arguments)

    return program_arguments


if __name__ == "__main__":
    main()

