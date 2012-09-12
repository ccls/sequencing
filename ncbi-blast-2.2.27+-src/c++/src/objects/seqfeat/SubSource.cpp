/* $Id: SubSource.cpp 371868 2012-08-13 15:10:25Z rafanovi $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  .......
 *
 * File Description:
 *   .......
 *
 * Remark:
 *   This code was originally generated by application DATATOOL
 *   using the following specifications:
 *   'seqfeat.asn'.
 */

// standard includes
#include <ncbi_pch.hpp>
#include <serial/enumvalues.hpp>

// generated includes
#include <objects/seqfeat/SubSource.hpp>

// generated classes

BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::

// destructor
CSubSource::~CSubSource(void)
{
}


void CSubSource::GetLabel(string* str) const
{
    typedef SStaticPair<TSubtype, const char*> TPair;
    static const TPair sc_Pairs[] = {
        { eSubtype_chromosome,              "chromosome" },
        { eSubtype_map,                     "map" },
        { eSubtype_clone,                   "clone" },
        { eSubtype_subclone,                "subclone" },
        { eSubtype_haplotype,               "haplotype" },
        { eSubtype_genotype,                "genotype" },
        { eSubtype_sex,                     "sex" },
        { eSubtype_cell_line,               "cell_line" },
        { eSubtype_cell_type,               "cell_type" },
        { eSubtype_tissue_type,             "tissue_type" },
        { eSubtype_clone_lib,               "clone_lib" },
        { eSubtype_dev_stage,               "dev_stage" },
        { eSubtype_frequency,               "frequency" },
        { eSubtype_germline,                "germline" },
        { eSubtype_rearranged,              "rearranged" },
        { eSubtype_lab_host,                "lab_host" },
        { eSubtype_pop_variant,             "pop_variant" },
        { eSubtype_tissue_lib,              "tissue_lib" },
        { eSubtype_plasmid_name,            "plasmid_name" },
        { eSubtype_transposon_name,         "transposon_name" },
        { eSubtype_insertion_seq_name,      "insertion_seq_name" },
        { eSubtype_plastid_name,            "plastid_name" },
        { eSubtype_country,                 "country" },
        { eSubtype_segment,                 "segment" },
        { eSubtype_endogenous_virus_name,   "endogenous_virus_name" },
        { eSubtype_transgenic,              "transgenic" },
        { eSubtype_environmental_sample,    "environmental_sample" },
        { eSubtype_isolation_source,        "isolation_source" },
        { eSubtype_lat_lon,                 "lat_lon" },
        { eSubtype_collection_date,         "collection_date" },
        { eSubtype_collected_by,            "collected_by" },
        { eSubtype_identified_by,           "identified_by" },
        { eSubtype_fwd_primer_seq,          "fwd_primer_seq" },
        { eSubtype_rev_primer_seq,          "rev_primer_seq" },
        { eSubtype_fwd_primer_name,         "fwd_primer_name" },
        { eSubtype_rev_primer_name,         "rev_primer_name" },
        { eSubtype_other,                   "other" }
    };
    typedef CStaticArrayMap<TSubtype, const char*> TSubtypeMap;
    DEFINE_STATIC_ARRAY_MAP(TSubtypeMap, sc_Map, sc_Pairs);

    *str += '/';
    TSubtypeMap::const_iterator iter = sc_Map.find(GetSubtype());
    if (iter != sc_Map.end()) {
        *str += iter->second;
    } else {
        *str += "unknown";
    }
    *str += '=';
    *str += GetName();
    if (IsSetAttrib()) {
        *str += " (";
        *str += GetAttrib();
        *str += ")";
    }
}


CSubSource::TSubtype CSubSource::GetSubtypeValue(const string& str)
{
    string name = NStr::TruncateSpaces(str);
    NStr::ToLower(name);
    replace(name.begin(), name.end(), '_', '-');

    return ENUM_METHOD_NAME(ESubtype)()->FindValue(str);
}


string CSubSource::GetSubtypeName(CSubSource::TSubtype stype)
{
    if (stype == CSubSource::eSubtype_other) {
        return "note";
    } else {
        return ENUM_METHOD_NAME(ESubtype)()->FindName(stype, true);
    }
}




bool CSubSource::NeedsNoText(const TSubtype& subtype)
{
    if (subtype == eSubtype_germline
        || subtype == eSubtype_rearranged
        || subtype == eSubtype_transgenic
        || subtype == eSubtype_environmental_sample
        || subtype == eSubtype_metagenomic) {
        return true;
    } else {
        return false;
    }
}


const string sm_LegalMonths [] = {
  "Jan",
  "Feb",
  "Mar",
  "Apr",
  "May",
  "Jun",
  "Jul",
  "Aug",
  "Sep",
  "Oct",
  "Nov",
  "Dec",
};


const int sm_daysPerMonth [] = {
  31,
  28,
  31,
  30,
  31,
  30,
  31,
  31,
  30,
  31,
  30,
  31
};

CRef<CDate> CSubSource::DateFromCollectionDate (const string& str) THROWS((CException))
{
    if (NStr::IsBlank(str)) {
        NCBI_THROW (CException, eUnknown,
                        "collection-date string is blank");
    }
    size_t pos = NStr::Find(str, "-");
    string year = "";
    string month = "";
    string day = "";
    int dpm = 0;

    if (pos == string::npos) {
        year = str;
    } else {
        size_t pos2 = NStr::Find(str, "-", pos + 1);
        if (pos2 == string::npos) {
            month = str.substr(0, pos);
            year = str.substr(pos + 1);
            if (NStr::IsBlank(month)) {
                NCBI_THROW (CException, eUnknown,
                                "collection-date string is improperly formatted");
            }
        } else {
            day = str.substr(0, pos);
            month = str.substr(pos + 1, pos2 - pos - 1);
            year = str.substr(pos2 + 1);
            if (NStr::IsBlank(month) || NStr::IsBlank(day)) {
                NCBI_THROW (CException, eUnknown,
                                "collection-date string is improperly formatted");
            }
        }
    }

    int day_val = 0;
    if (!NStr::IsBlank(day)) {
        try {
            day_val = NStr::StringToInt (day);
            if (day_val < 1 || day_val > 31) {
                NCBI_THROW (CException, eUnknown,
                                "collection-date string has invalid day value");
            }
        } catch ( const exception& ) {
            // threw exception while converting to int
            NCBI_THROW (CException, eUnknown,
                            "collection-date string is improperly formatted");
        }
    }

    int month_val = 0;
    if (!NStr::IsBlank(month)) {
        for (size_t i = 0; i < sizeof(sm_LegalMonths) / sizeof(string); i++) {
            if (NStr::Equal(month, sm_LegalMonths[i])) {
                dpm = sm_daysPerMonth [i];
                month_val = i + 1;
                break;
            }
        }
        if (month_val == 0) {
            NCBI_THROW (CException, eUnknown,
                            "collection-date string has invalid month");
        }
    }

    if (NStr::IsBlank(year)) {
        NCBI_THROW (CException, eUnknown,
                        "collection-date string is improperly formatted");
    } 

    int year_val = 0;
    try {
        year_val = NStr::StringToInt (year);
    } catch ( const exception& ) {
        // threw exception while converting to int
        NCBI_THROW (CException, eUnknown,
                        "collection-date string is improperly formatted");
    }

    if (year_val < 1700 || year_val >= 2100) {
        NCBI_THROW (CException, eUnknown,
                        "collection-date year is out of range");
    }

    if (day_val > 0 && dpm > 0 && day_val > dpm) {
        if (month_val != 2 || day_val != 29 || (year_val % 4) != 0) {
            NCBI_THROW (CException, eUnknown,
                            "collection-date day is greater than monthly maximum");
        }
    }

    CRef<CDate> date(new CDate);

    date->SetStd().SetYear (year_val);
    if (month_val > 0) {
        date->SetStd().SetMonth (month_val);
    }
    if (day_val > 0) {
        date->SetStd().SetDay (day_val);
    }
    
    return date;
}


// =============================================================================
//                                 Country Names
// =============================================================================


// legal country names
const string CCountries::sm_Countries[] = {
    "Afghanistan",
    "Albania",
    "Algeria",
    "American Samoa",
    "Andorra",
    "Angola",
    "Anguilla",
    "Antarctica",
    "Antigua and Barbuda",
    "Arctic Ocean",
    "Argentina",
    "Armenia",
    "Aruba",
    "Ashmore and Cartier Islands",
    "Atlantic Ocean",
    "Australia",
    "Austria",
    "Azerbaijan",
    "Bahamas",
    "Bahrain",
    "Baker Island",
    "Baltic Sea",
    "Bangladesh",
    "Barbados",
    "Bassas da India",
    "Belarus",
    "Belgium",
    "Belize",
    "Benin",
    "Bermuda",
    "Bhutan",
    "Bolivia",
    "Borneo",
    "Bosnia and Herzegovina",
    "Botswana",
    "Bouvet Island",
    "Brazil",
    "British Virgin Islands",
    "Brunei",
    "Bulgaria",
    "Burkina Faso",
    "Burundi",
    "Cambodia",
    "Cameroon",
    "Canada",
    "Cape Verde",
    "Cayman Islands",
    "Central African Republic",
    "Chad",
    "Chile",
    "China",
    "Christmas Island",
    "Clipperton Island",
    "Cocos Islands",
    "Colombia",
    "Comoros",
    "Cook Islands",
    "Coral Sea Islands",
    "Costa Rica",
    "Cote d'Ivoire",
    "Croatia",
    "Cuba",
    "Curacao",
    "Cyprus",
    "Czech Republic",
    "Democratic Republic of the Congo",
    "Denmark",
    "Djibouti",
    "Dominica",
    "Dominican Republic",
    "East Timor",
    "Ecuador",
    "Egypt",
    "El Salvador",
    "Equatorial Guinea",
    "Eritrea",
    "Estonia",
    "Ethiopia",
    "Europa Island",
    "Falkland Islands (Islas Malvinas)",
    "Faroe Islands",
    "Fiji",
    "Finland",
    "France",
    "French Guiana",
    "French Polynesia",
    "French Southern and Antarctic Lands",
    "Gabon",
    "Gambia",
    "Gaza Strip",
    "Georgia",
    "Germany",
    "Ghana",
    "Gibraltar",
    "Glorioso Islands",
    "Greece",
    "Greenland",
    "Grenada",
    "Guadeloupe",
    "Guam",
    "Guatemala",
    "Guernsey",
    "Guinea-Bissau",
    "Guinea",
    "Guyana",
    "Haiti",
    "Heard Island and McDonald Islands",
    "Honduras",
    "Hong Kong",
    "Howland Island",
    "Hungary",
    "Iceland",
    "India",
    "Indian Ocean",
    "Indonesia",
    "Iran",
    "Iraq",
    "Ireland",
    "Isle of Man",
    "Israel",
    "Italy",
    "Jamaica",
    "Jan Mayen",
    "Japan",
    "Jarvis Island",
    "Jersey",
    "Johnston Atoll",
    "Jordan",
    "Juan de Nova Island",
    "Kazakhstan",
    "Kenya",
    "Kerguelen Archipelago",
    "Kingman Reef",
    "Kiribati",
    "Kosovo",
    "Kuwait",
    "Kyrgyzstan",
    "Laos",
    "Latvia",
    "Lebanon",
    "Lesotho",
    "Liberia",
    "Libya",
    "Liechtenstein",
    "Lithuania",
    "Luxembourg",
    "Macau",
    "Macedonia",
    "Madagascar",
    "Malawi",
    "Malaysia",
    "Maldives",
    "Mali",
    "Malta",
    "Marshall Islands",
    "Martinique",
    "Mauritania",
    "Mauritius",
    "Mayotte",
    "Mediterranean Sea",
    "Mexico",
    "Micronesia",
    "Midway Islands",
    "Moldova",
    "Monaco",
    "Mongolia",
    "Montenegro",
    "Montserrat",
    "Morocco",
    "Mozambique",
    "Myanmar",
    "Namibia",
    "Nauru",
    "Navassa Island",
    "Nepal",
    "Netherlands",
    "New Caledonia",
    "New Zealand",
    "Nicaragua",
    "Niger",
    "Nigeria",
    "Niue",
    "Norfolk Island",
    "North Korea",
    "North Sea",
    "Northern Mariana Islands",
    "Norway",
    "Oman",
    "Pacific Ocean",
    "Pakistan",
    "Palau",
    "Palmyra Atoll",
    "Panama",
    "Papua New Guinea",
    "Paracel Islands",
    "Paraguay",
    "Peru",
    "Philippines",
    "Pitcairn Islands",
    "Poland",
    "Portugal",
    "Puerto Rico",
    "Qatar",
    "Republic of the Congo",
    "Reunion",
    "Romania",
    "Ross Sea",
    "Russia",
    "Rwanda",
    "Saint Helena",
    "Saint Kitts and Nevis",
    "Saint Lucia",
    "Saint Pierre and Miquelon",
    "Saint Vincent and the Grenadines",
    "Samoa",
    "San Marino",
    "Sao Tome and Principe",
    "Saudi Arabia",
    "Senegal",
    "Serbia",
    "Seychelles",
    "Sierra Leone",
    "Singapore",
    "Sint Maarten",
    "Slovakia",
    "Slovenia",
    "Solomon Islands",
    "Somalia",
    "South Africa",
    "South Georgia and the South Sandwich Islands",
    "South Korea",
    "Southern Ocean",
    "Spain",
    "Spratly Islands",
    "Sri Lanka",
    "Sudan",
    "Suriname",
    "Svalbard",
    "Swaziland",
    "Sweden",
    "Switzerland",
    "Syria",
    "Taiwan",
    "Tajikistan",
    "Tanzania",
    "Tasman Sea",
    "Thailand",
    "Togo",
    "Tokelau",
    "Tonga",
    "Trinidad and Tobago",
    "Tromelin Island",
    "Tunisia",
    "Turkey",
    "Turkmenistan",
    "Turks and Caicos Islands",
    "Tuvalu",
    "Uganda",
    "Ukraine",
    "United Arab Emirates",
    "United Kingdom",
    "Uruguay",
    "USA",
    "Uzbekistan",
    "Vanuatu",
    "Venezuela",
    "Viet Nam",
    "Virgin Islands",
    "Wake Island",
    "Wallis and Futuna",
    "West Bank",
    "Western Sahara",
    "Yemen",
    "Zambia",
    "Zimbabwe"
};

const string CCountries::sm_Former_Countries[] = {
    "Belgian Congo",
    "British Guiana",
    "Burma",
    "Czechoslovakia",
    "Netherlands Antilles",
    "Korea",
    "Serbia and Montenegro",
    "Siam",
    "USSR",
    "Yugoslavia",
    "Zaire"
};

bool CCountries::IsValid(const string& country)
{
    string name = country;
    size_t pos = country.find(':');

    if ( pos != string::npos ) {
        name = country.substr(0, pos);
    }

    // try current countries
    const string *begin = sm_Countries;
    const string *end = &(sm_Countries[sizeof(sm_Countries) / sizeof(string)]);

    return find(begin, end, name) != end;
}


bool CCountries::IsValid(const string& country, bool& is_miscapitalized)
{
    string name = country;
    size_t pos = country.find(':');

    if ( pos != string::npos ) {
        name = country.substr(0, pos);
    }

    is_miscapitalized = false;
    // try current countries
    size_t num_countries = sizeof(sm_Countries) / sizeof(string);
    for (size_t i = 0; i < num_countries; i++) {
        if (NStr::EqualNocase (name, sm_Countries[i])) {
            if (!NStr::EqualCase(name, sm_Countries[i])) {
                is_miscapitalized = true;
            }
            return true;
        }
    }
    return false;
}


bool CCountries::WasValid(const string& country)
{
    string name = country;
    size_t pos = country.find(':');

    if ( pos != string::npos ) {
        name = country.substr(0, pos);
    }

    // try formerly-valid countries
    const string *begin = sm_Former_Countries;
    const string *end = &(sm_Former_Countries[sizeof(sm_Former_Countries) / sizeof(string)]);

    return find(begin, end, name) != end;
}


bool CCountries::WasValid(const string& country, bool& is_miscapitalized)
{
    string name = country;
    size_t pos = country.find(':');

    if ( pos != string::npos ) {
        name = country.substr(0, pos);
    }

    is_miscapitalized = false;
    // try formerly-valid countries
    size_t num_countries = sizeof(sm_Former_Countries) / sizeof(string);
    for (size_t i = 0; i < num_countries; i++) {
        if (NStr::EqualNocase (name, sm_Former_Countries[i])) {
            if (!NStr::EqualCase(name, sm_Former_Countries[i])) {
                is_miscapitalized = true;
            }
            return true;
        }
    }
    return false;
}


END_objects_SCOPE // namespace ncbi::objects::

END_NCBI_SCOPE

/* Original file checksum: lines: 65, chars: 1891, CRC32: 7724f0c5 */