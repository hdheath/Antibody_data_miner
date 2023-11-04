# This function will load the Neuromab database as df

def NeuromabDF():
    
    # Grab the dataframe from the Neuromab site
    data = requests.get("https://neuromab.ucdavis.edu/catalog-download.cfm").content
    df = pd.read_csv(BytesIO(data))
    
    return df

# This is the main call function that will loop each pdf name of the Neuromab df, and also calls the Accession Search and Peptide Slice
def PDFminer():
    
    df = NeuromabDF()
    
    df['Data Sheet File Name'] = 0
    df['NeuroMab Clone'] = 0
    df["Accession Number"] = 0
    df['Amino Acid Range'] = 0
    df['Organism'] = 0
    df['Full Amino Seq'] = 0
    df['Peptide Sequence'] = 0
    
    
    # List of Antibodies lacking information (will make into dict.)
    NA_list = ['N286_74.pdf', 'N175_9.pdf', 'A12_18.pdf', 'N358_68.pdf', 
               'N297_59.pdf', 'N257_25.pdf', 'N259_48.pdf','N262_16.pdf', 
               'N144_14.pdf','N159_5.pdf']
    
    # Iterate over all pdf files in the directory
    for i, row in df.iterrows(): 
    
        pdf_name = row['DataSheetFileName']
        pdf_name = str(pdf_name)
        
        accession_number = row['AccessionNum']
        accession_number = str(accession_number).strip()
        
        
        if pdf_name not in NA_list:
            
            text = PDFConverter(pdf_name,i)
            
            df.at[i, 'Data Sheet File Name'] = (pdf_name)
            df.at[i, 'NeuroMab Clone'] = NeuromabScan(text,accession_number,pdf_name)
            df.at[i, "Accession Number"] = AccessionScan(text,accession_number,pdf_name)
            df.at[i,  'Amino Acid Range'] = AminoScan(text,accession_number,pdf_name)
            aa_range = df.at[i,  'Amino Acid Range']
            df.at[i,  'Organism'] = OrganismScan(text,accession_number,pdf_name)
            df.at[i, 'Full Amino Seq'] = AccessionSearch(i,row,accession_number,df,pdf_name)
            df.at[i, 'Peptide Sequence'] = SplicePeptide(i,row,aa_range,df,accession_number,pdf_name)

        else:
            df.at[i, 'Data Sheet File Name'] = (pdf_name)
            df.at[i, 'NeuroMab Clone'] = 'NA'
            df.at[i, "Accession Number"] = 'NA'
            df.at[i, 'Amino Acid Range'] = 'NA'
            df.at[i, 'Organism'] = 'NA'
            df.at[i, 'Full Amino Seq'] = 'NA'
            df.at[i, 'Peptide Sequence'] = 'NA'
            
            print ('PDFMiner not executed on' + ' ' + pdf_name + ' ' + accession_number)
            
    return df

# This function converts a url into text

def PDFConverter(
    pdf_name,i: str,
)-> str:
    
    # Set the beginning url name to access indivdual pdfs
    path = 'https://neuromab.ucdavis.edu/datasheet/'

    # Get the file path
    pdf_path = path + pdf_name

    # Download the PDF file
    response = requests.get(pdf_path)

    # Save the PDF file to a local file
    with open(pdf_name, 'wb') as file:
        file.write(response.content)

    # Open the PDF file in binary mode
    with open(pdf_name, 'rb') as file:
        
        try:

            # Create a PDF parser object
            parser = PDFParser(file)

            # Create a PDF document object
            doc = PDFDocument(parser)

            # Create a PDF resource manager object
            resource_manager = PDFResourceManager()

            # Create a buffer for the extracted text
            output = StringIO()

            # Create a text converter object
            converter = TextConverter(
                resource_manager, output, laparams=LAParams()
            )

            # Create a PDF interpreter object
            interpreter = PDFPageInterpreter(resource_manager, converter)

            # Process each page contained in the pdf file
            for page in PDFPage.get_pages(file, caching=True, check_extractable=True):
                interpreter.process_page(page)

            # Close the converter
            converter.close()

            # Get the text from the buffer
            text = output.getvalue()

            # Close the buffer
            output.close()
            
        except:
            print('PDFconverter failed for' + (' ') + pdf_name)
            text = 'NA'
            
    return text

# This function indexes a string for the object after "NeuroMab Clone" - which should be the Neuormab ID (used to cross-check with the Neuromab database)
def NeuromabScan(
    text,accession_number,pdf_name: str,
) -> str:

    # Search for the string "NeuroMab clone" in the text
    if "NeuroMab clone" in text:
        
        # Get the index of the first occurrence of the string
        index = text.index("NeuroMab clone")

        # Get the string that comes immediately after the string "NeuroMab clone"
        NeuroMab = text[index + len("NeuroMab clone"):].split()[0]
            
        # Remove any whitespace or parenthesis from the value
        NeuroMab = NeuroMab.strip().strip("()")
        return NeuroMab

        
    else:
        NeuroMab = "----"
        print('NeuroMab failed for:' + ' ' + pdf_name + ' ' + accession_number)
        return NeuroMab


# This function indexes a string for the organism type from the NeuroMab pdfs
def OrganismScan(
    text,accession_number,pdf_name: str,
) -> str:

    # Search for the organism string in the text
    if ("of rat" in text) or ('Mouse=' and "Human:" in text) or ('Mouse:' and "Human:" in text):
        return "Rat"
    
    elif ("of human" in text) or ('Mouse:' and "Rat:" in text) or ('Mouse=' and "Rat:" in text):
        return "Human"
    
    elif ("of mouse" in text) or ("Rat:" and "Human:" in text):
        return "Mouse"
        
    elif ("of goldfish" in text):      
        return "Goldfish"
        
    elif "of jellyfish" in text:      
        return "Jellyfish"
    
    elif "of zebrafish" in text:       
        return "Zebrafish"
        
    else:
        print('OrganismScan failed for:' + ' ' + pdf_name + ' ' + accession_number)
        return '----'


# This function cleans the string of the Amino Range found in the pdf text
def AminoStringCleaner(
    Amino: str,
) -> str:
    # Remove any whitespace or parenthesis from the value
        Amino = ('(') + Amino.strip().strip("()").strip(",") + (')')
        return Amino

# This function looks for the string after 'protein' for the amino acid range if the search for 'acids' in the text provided the wrong information
def ProteinScan(
    text: str,
) -> str:
    # Get the index of the first occurrence of the string
    index = text.index("protein")

    # Get the string that comes immediately after the string "acids"
    Amino = text[index + len("protein"):].split()[0]
    
    # Clean the string 
    Amino = AminoStringCleaner(Amino)
    
    # Correct if range not found 
    if ('-') not in Amino:
        Amino = 'NA'
    
    return Amino 

# This function looks for the string after 'peptide' for the amino acid range if the search for 'acids' in the text provided the wrong information
def PeptideScan(
    text: str,
) -> str:
    # Get the index of the first occurrence of the string
    index = text.index("peptide")

    # Get the string that comes immediately after the string "acids"
    Amino = text[index + len("peptide"):].split()[0]
    
    # Clean the string 
    Amino = AminoStringCleaner(Amino)
    
    # Correct if range not found 
    if ('-') not in Amino:
        Amino == 'NA'
        
    
    return Amino

# This function indexes a string for the object after "acids" for the amino acid range applicable to the Neuromab Clone from a whole peptide sequence

def AminoScan(
    text,accession_number,pdf_name: str,
) -> str:
    
    Amino = ''
    
    if "1D8" in pdf_name:
        return '(904-917)'
    
    elif 'acids' in text:
        # Get the index of the first occurrence of the string
        index = text.index("amino acids")

        # Get the string that comes immediately after the string "acids"
        Amino = text[index + len("amino acids"):].split()[0]
        
        # Clean the string 
        Amino = AminoStringCleaner(Amino)
        
        # Special case for 1D8
        if Amino == ('(Leucine-rich)'):
            Amino = '(904-917)'
        
        # if the Amino acid range is not found correctly from 'acids', 
        #look for protein or peptide
        elif ('-' not in Amino) or (len(Amino) < 4):
            if 'protein' in text:
                Amino = ProteinScan(text)
            
            elif 'peptide' in text:
                Amino = PeptideScan(text)
    
    # In the event that 'acids' is not in the text, search for 
    # protein of peptide instead 
    elif ('-' not in Amino) or (len(Amino) < 4):
        if 'protein' in text:
                Amino = ProteinScan(text)
            
        elif 'peptide' in text:
            Amino = PeptideScan(text)
            
    else:
        print('AminoScan failed for:' + ' ' + pdf_name + ' ' + accession_number)
        return 'NA'
    
    return Amino


# This function searches for the accession number from the pdfs (to cross check with Neuromab database/fill holes)
def AccessionScan(
    text,accession_number,pdf_name: str,
) -> str:

    # Search for the string "accession number" in the text
    if "number" in text:
        # Get the index of the first occurrence of the string
        index = text.index("number")

        # Get the string that comes immediately after the string "accession number"
        AN = text[index + len("number"):].split()[0]
            
        # Remove any whitespace or parenthesis from the value
        AN = AN.strip().strip("()")
        
        if (AN == 'P16389),'):
            AN = 'P16389'
        
        elif (AN == 'O95180),'):
            AN = 'O95180'
            
    else:
        print('AccessionScan failed for:' + ' ' + pdf_name + ' ' + accession_number)
        AN = "-----"
    
    return AN


# This function searches the Uniprot database using a Uniprot ID

def UniprotIDSearch(
    path,accession_number,json: str
) -> str:
    
    # Combine the shared path and the URL information to create the full URL
    url = path + accession_number + json
    
    # Set a key to search for 
    key = 'sequence'
    
    # Send a request to the URL and retrieve the HTML content

    try: 
        json_file = (requests.get(url).json())
        
        if key in json_file:
            # If there is a response and key in HTML grab the sequence
            return (json_file[key])
        else:
            return 'NA'
            
    except:
        print('UniprotIDSearch failed for:' + ' ' + url)
        return 'NA'    

# This function uses an NCBI ID to search the web for the correct accession
def NCBIiDSearch(
    path,accession_number,json: str
) -> str:
    
    # set url parameters 
    path2 = 'https://rest.uniprot.org/uniprotkb/stream?compressed=false&download=false&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&format=tsv&query=%28'
    end = '%29'
    url2 = (path2) + (accession_number) + (end)
    try:
        # Send a request to the URL and retrieve the Uniprot Accession
        url_info2 = requests.get(url2)
        url_info2 = url_info2.text.strip('').split('\n')[1][0:6]
        accession_number=url_info2

        # Send the new accession to the Uniprot search function 
        sequence = UniprotIDSearch(path,accession_number,json)
    except:
        print('NCBIidSearch failed for:' + ' ' + url2)
        sequence = 'NA'
    
    return sequence


# This function takes the accession numbers from the Neuromab database, and calls upon the correct function to search the web for the peptide sequence
def AccessionSearch(
    i,row,accession_number,df,pdf_name : str
) -> str:
    
    path = 'https://rest.uniprot.org/uniprotkb/'
    json = '.json'
    
    isotype_list = ['P18507-2', 'Q92913-2', 'Q9D8C3-2','Q8BR86-3','P22462-3','P10827-2',
                    'Q96PU8-2','Q96PU8-5','Q96PU8-8','Q8NFF2-3','P10827-2','Q8NFF2-3',
                    'Q63563-2','Q8C437-3','Q14722-3','P04896-2','P35439-2','Q92913-2',
                    'Q9D8C3-2', 'Q8C437-3','Q62634', 'Q16635-3']
    
    retry_list = ['P10827-2','Q63563-2']
        
    # if an Uniprot ID is found, search using 
    #UniprotIDSearch()
    if (len(accession_number) < 8) and (len(accession_number) > 3): 
        return UniprotIDSearch(path,accession_number,json)
                  
    # if a NCBI is found use another method to search 
    #for the sequence
    elif ("NP" in accession_number) or ("AA" in accession_number) or ("AB" in accession_number):
        return NCBIiDSearch(path,accession_number,json)
    
    elif (accession_number in isotype_list):
        return UniprotIDSearch(path,accession_number,json)
    
    elif 'XM_342497' in accession_number:
        return ('{' + 'value' + ':' + 'MKQESAAQSTPPPSLSPAPSAQPSWEDGDPQALWIFGYGSLVWKPDFAYSDSRVGFVRGYSRRFWQGDTFHRGSDKMPGRVVTLLEDHEGCTWGVAYQVRGEQVSEALKYLNVREAVLGGYDTKEVTFYPQDTPDQPLTALAYVATPQNPGYLGPAPEEVIATQILACRGFSGHNLEYLLRLADFMQLCGPQAQDEHLEAIVDAVGSLLPCSYLSEQPLALI')
    elif 'XP_989315' in accession_number:
        return ('{' + 'value' + ':' + '“MMGRRGLPGHRPRPAQPPAPSPALPFPARQAQQHSPRHPAHAPQRAARRHYSAGQQKTSIGDFRLTTGSHGYHGHNSGGEHPIQTSPGLMQPFSIPVQITLQGGRRRQGRTALPASGKINGDPLKVHPKLPSSAGEDRAMLLGVAMMASSVLMFFLLGTTVLKPFMLRSPREESNCTTVHTHIVDDGLDFSFTCEGSCQDHGRSPCLQVFVNLSHSGQKVLLHYDDEAIRTNPKCFYTPKCHGDRDDLLNSVLDIKEFFDHNNGTFPCFYSPDGPLGVVLRKSGHKVVFHCLFWPLLTLLGGALIVGLVRLTQHLSFQCEQYSTVVRA"')
    
    # if the accession row is empty, use the acession
    # found from the pdfs
    elif (accession_number) or (len(accession_number) > 2) :
        try: 
            accession_number = df.at[i, 'Accession Number']
            
            # if an Uniprot ID is found, search using 
            #UniprotIDSearch()
            if (len(accession_number) < 8) and (len(accession_number) > 3):
                try:
                    # Fill in blank space with accession
                    df.at[i, 'AccessionNum'] == accession_number
                except:
                    print('Could not change df')
                
                return UniprotIDSearch(path,accession_number,json)

            # if a NCBI is found use another method to search 
            #for the sequence
            elif ("NP" in accession_number) or ("AA" in accession_number) or ("AB" in accession_number):     
                try:
                    # Fill in blank space with accession
                    df.at[i, 'AccessionNum'] = accession_number
                
                except:
                    print('Could not change df')
                
                return NCBIiDSearch(path,accession_number,json)
        except:
            print('AccessionSearch failed for:' + ' ' + pdf_name + ' ' + accession_number) 
            return ("NA")
    elif (accession_number) in retry_list:
        accession_number = df.at[i, 'Accession Number']
        return UniprotIDSearch(path,accession_number,json)
    else:
        print('AccessionSearch failed for:' + ' ' + pdf_name + ' ' + accession_number) 
        return ("NA")    

# This function is for Neuromabs that did not have their amino range in the pdf, instead the spliced sequences are in the pdf themself (N106)

def SequenceInPDF(
    pdf_name,i: str,
) -> str:
    text = PDFConverter(pdf_name,i)
    
    # Look for sequence in the pdf
    if "sequence:" in text:  
        start = text.index("sequence:")
        end = text.index("EKKAH") + len("sequence:")
        substring = text[start:end]
        substring = substring[10:-2]
        substring = str(substring)
        return substring
        
    elif 'AB_10671176' in text:
        sequence = 'DENYHNERARKSRNRLS'
    
    elif 'NP_001020928' in text:
        sequence = 'KTTLYAFDELDEFPETSV'

    else:
        sequence = 'NA'
    return sequence

# This function adds the spliced peptide sequences to the dataframe
def SplicePeptide(
    i,row,aa_range,df,accession_number,pdf_name : str
) -> str:

    
    # Get the sequence string and pdf name 

    string = df.at[i, 'Full Amino Seq']
    data_file_name = df.at[i, 'DataSheetFileName']
    data_file_name = str(data_file_name)

    # Remove the '{'value':' substring from the string
    string = str(string)
    string = string[9:]
    
        # Look for sequence within pdf 
    if ('N106' in data_file_name) or ('N10_7' in data_file_name) or ('N15_' in pdf_name) or ('N116_14' in data_file_name):
        return SequenceInPDF(data_file_name,i)

    elif ('-' in aa_range):
        
        try: 
            
            # Split the range into start and end indices
            # Adjust for 0 based start 
            aa_range = aa_range.strip('()')
            start, end = aa_range.split("-")
            start = int(start)
            start = start +1
            end = int(end)
            end = end + 2


            # Splice the string and append it to the list
            return string[start:end]
        
        except:
            print ('SlicePeptide failed for' + ' ' + pdf_name + ' ' + accession_number) 
            return 'NA'
                                                                                                                 
    else:        
        print ('SlicePeptide failed for' + ' ' + pdf_name + ' ' + accession_number) 
        return "NA"

  
