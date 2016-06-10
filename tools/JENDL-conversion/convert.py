import endfgeant4converter
reload(endfgeant4converter)
egc = endfgeant4converter.ENDFGeant4Converter(inputdir = "JENDL-AN-2005-linux-patched")

egc.convertF0()
egc.convert()
