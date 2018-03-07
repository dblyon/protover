aa_reference = {"A": {"name": "Alanine",
                      "tlc": "Ala",
                      "comp": {"C": 3,
                               "H": 7,
                               "N": 1,
                               "O": 2,
                               "S": 0}                 
                      },
                "R": {"name": "Arginine",
                      "tlc": "Arg",
                      "comp": {"C": 6,
                               "H": 14,
                               "N": 4,
                               "O": 2,
                               "S": 0}
                      },
                "N": {"name": "Asparagine",
                      "tlc": "Asn",
                      "comp": {"C": 4,
                               "H": 8,
                               "N": 2,
                               "O": 3,
                               "S": 0}
                      },
                "D": {"name": "Asparticacid",
                      "tlc": "Asp",
                      "comp": {"C": 4,
                               "H": 7,
                               "N": 1,
                               "O": 4,
                               "S": 0}
                      },          
                "C": {"name": "Cysteine",
                      "tlc": "Cys",
                      "comp": {"C": 3,
                               "H": 7,
                               "N": 1,
                               "O": 2,
                               "S": 1}
                      },
                "G": {"name": "Gylcine",
                      "tlc": "Gly",
                      "comp": {"C": 2,
                               "H": 5,
                               "N": 1,
                               "O": 2,
                               "S": 0}
                      },
                "Q": {"name": "Glutamine",
                      "tlc": "Gln",
                      "comp": {"C": 5,
                               "H": 10,
                               "N": 2,
                               "O": 3,
                               "S": 0}
                      },
                "E": {"name": "Glutamicacid",
                      "tlc": "Glu",
                      "comp": {"C": 5,
                               "H": 9,
                               "N": 1,
                               "O": 4,
                               "S": 0}
                      },
                "H": {"name": "Histidine",
                      "tlc": "His",
                      "comp": {"C": 6,
                               "H": 9,
                               "N": 3,
                               "O": 2,
                               "S": 0}
                      },
                "I": {"name": "Isoleucine",
                      "tlc": "Ile",
                      "comp": {"C": 6,
                               "H": 13,
                               "N": 1,
                               "O": 2,
                               "S": 0}
                      },
                "L": {"name": "Leucine",
                      "tlc": "Leu",
                      "comp": {"C": 6,
                               "H": 13,
                               "N": 1,
                               "O": 2,
                               "S": 0}
                      },
                "K": {"name": "Lysine",
                      "tlc": "Lys",
                      "comp": {"C": 6,
                               "H": 14,
                               "N": 2,
                               "O": 2,
                               "S": 0}
                      },
                "M": {"name": "Methionine",
                      "tlc": "Met",
                      "comp": {"C": 5,
                               "H": 11,
                               "N": 1,
                               "O": 2,
                               "S": 1}
                      },
                "F": {"name": "Phenylalanine",
                      "tlc": "Phe",
                      "comp": {"C": 9,
                               "H": 11,
                               "N": 1,
                               "O": 2,
                               "S": 0}
                      },
                "P": {"name": "Proline",
                      "tlc": "Pro",
                      "comp": {"C": 5,
                               "H": 9,
                               "N": 1,
                               "O": 2,
                               "S": 0}
                      },
                "S": {"name": "Serine",
                      "tlc": "Ser",
                      "comp": {"C": 3,
                               "H": 7,
                               "N": 1,
                               "O": 3,
                               "S": 0}
                      },
                "T": {"name": "Threonine",
                      "tlc": "Thr",
                      "comp": {"C": 4,
                               "H": 9,
                               "N": 1,
                               "O": 3,
                               "S": 0}
                      },
                "W": {"name": "Tryptophan",
                      "tlc": "Trp",
                      "comp": {"C": 11,
                               "H": 12,
                               "N": 2,
                               "O": 2,
                               "S": 0}
                      },
                "Y": {"name": "Tyrosine",
                      "tlc": "Tyr",
                      "comp": {"C": 9,
                               "H": 11,
                               "N": 1,
                               "O": 3,
                               "S": 0}
                      },
                "V": {"name": "Valine",
                      "tlc": "Val",
                      "comp": {"C": 5,
                               "H": 11,
                               "N": 1,
                               "O": 2,
                               "S": 0}
                      },
                "U": {"name": "Selenocysteine",
                      "tlc": "Sec",
                      "comp": {"C": 3,
                               "H": 7,
                               "N": 1,
                               "O": 2,
                               "S": 0,
                               "Se": 1}
                      },                      
                "Mx": {"name": "Methionine_oxidated",
                       "tlc": "Metox",
                       "comp": {"C": 5,
                               "H": 11,
                               "N": 1,
                               "O": 3,
                               "S": 1}
                       },
                "X": {"name": "any",                  # average of the chemical composition of the 20 AAs
                      "tlc": "any",
                      "comp": {"C": 5,
                               "H": 10,
                               "N": 1,
                               "O": 2,
                               "S": 0}
                       },
                "B": {"name": "Aspartate/Asparagine", # converts to Asparagine
                      "tlc": "Asp/Asn",
                      "comp": {"C": 4,
                               "H": 8,
                               "N": 2,
                               "O": 3,
                               "S": 0}
                      },       
                "Z": {"name": "Glutamate/Glutamine", # converts to Glutamine
                      "tlc": "Asp/Asn",
                      "comp": {"C": 5,
                               "H": 10,
                               "N": 2,
                               "O": 3,
                               "S": 0}
                      } 
                }

def oneto3(code):
    """
    enter one letter amino acid code
    returns three letter code equivalent
    """
    code = code.upper().strip()
    return aa_reference[code]["tlc"]

def onetoname(code):
    """
    enter one letter amino acid code
    returns entire name equivalent
    """
    code = code.upper().strip()
    return aa_reference[code]["name"]

def threeto1(code):
   """
   enter three letter amino acid code
   returns one letter code equivalent
   """
   code = code.capitalize().strip()
   for key in aa_reference:
       if aa_reference[key]["tlc"] == code:
           return key
        


