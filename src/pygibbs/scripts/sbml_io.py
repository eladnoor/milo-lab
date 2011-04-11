import libsbml

filename = '../data/thermodynamics/yeast_4.02.xml'

reader = libsbml.SBMLReader()
document = reader.readSBMLFromFile(filename)
if document.getNumErrors():
    raise Exception('cannot read SBML model from file %s due to error: %s' % 
                    (filename, document.getError(0).getMessage()))

model = document.getModel()    
for i in xrange(model.getNumSpeciesTypes()):
    species_type = model.getSpeciesType(i)
    ann = species_type.getAnnotation()
    if ann:
        print i, ann.toXMLString()