import logging
import os
import tempfile
import shutil
from PyPDF2 import PdfFileWriter, PdfFileReader

logger=logging.getLogger("pdfedit")

def AddPdfKeys(pdf_file, keydict):
    """
    Add a dictionary of key-value pairs to a pdf file.
    Adds metadata to PDF files so you know what generated them.
    """
    multifiles=list()
    if isinstance(pdf_file, list):
        multifiles=pdf_file
    else:
        multifiles=[pdf_file]
    for of in multifiles:
        outstream=tempfile.NamedTemporaryFile(mode="wb", delete=False)
        logger.debug("AddPdfKeys created {0}".format(outstream.name))
        CopyPdfKeys(of, outstream, keydict)
        outstream.close()
        shutil.copyfile(outstream.name, of)
        os.remove(outstream.name)

def CopyPdfKeys(pdf_file, outstream, keydict):
    pdfDict={u"/{0}".format(k) : v for (k, v) in keydict.items()}
    infile=PdfFileReader(pdf_file)
    outfile=PdfFileWriter()
    outfile.addMetadata(pdfDict)
    for page in range(infile.getNumPages()):
        outfile.addPage(infile.getPage(page))
    outfile.write(outstream)
