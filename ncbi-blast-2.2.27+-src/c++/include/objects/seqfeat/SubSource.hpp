/* $Id: SubSource.hpp 371868 2012-08-13 15:10:25Z rafanovi $
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
 */

/// @SubSource.hpp
/// User-defined methods of the data storage class.
///
/// This file was originally generated by application DATATOOL
/// using the following specifications:
/// 'seqfeat.asn'.
///
/// New methods or data members can be added to it if needed.
/// See also: SubSource_.hpp


#ifndef OBJECTS_SEQFEAT_SUBSOURCE_HPP
#define OBJECTS_SEQFEAT_SUBSOURCE_HPP


// generated includes
#include <objects/seqfeat/SubSource_.hpp>

// generated classes


// other includes
#include <objects/general/Date.hpp>
#include <objects/general/Date_std.hpp>

BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::
class CDate;
class CDate_std;

/////////////////////////////////////////////////////////////////////////////
class NCBI_SEQFEAT_EXPORT CSubSource : public CSubSource_Base
{
    typedef CSubSource_Base Tparent;
public:
    // constructor
    CSubSource(void);
    CSubSource(TSubtype subtype, const TName& name);
    CSubSource(const string& subtype, const TName& name);

    // destructor
    ~CSubSource(void);

    void GetLabel(string* str) const;

    // convert subtype from string to enum.
    static TSubtype GetSubtypeValue(const string& str);

	// get name for subsource
    static string GetSubtypeName(CSubSource::TSubtype stype);

	// identify whether subsource value should be blank
	static bool NeedsNoText (const TSubtype& subtype);

	// read collection date from string
    static CRef<CDate> DateFromCollectionDate (const string& str) THROWS((CException));

private:
    // Prohibit copy constructor and assignment operator
    CSubSource(const CSubSource& value);
    CSubSource& operator=(const CSubSource& value);

};

/////////////////// CSubSource inline methods

// constructor
inline
CSubSource::CSubSource(void)
{
}

inline
CSubSource::CSubSource(TSubtype subtype, const TName& name)
{
    SetSubtype(subtype);
    SetName(name);
}

inline
CSubSource::CSubSource(const string& subtype, const TName& name)
{
    SetSubtype(GetSubtypeValue(subtype));
    SetName(name);
}


/////////////////// end of CSubSource inline methods


// =============================================================================
//                 Country Names (legal values found in country subtype)
// =============================================================================


class NCBI_SEQFEAT_EXPORT CCountries
{
public:
    static bool IsValid(const string& country);
    static bool IsValid(const string& country, bool& is_miscapitalized);
    static bool WasValid(const string& country);
    static bool WasValid(const string& country, bool& is_miscapitalized);

private:
    static const string sm_Countries[];
    static const string sm_Former_Countries[];
};




END_objects_SCOPE // namespace ncbi::objects::

END_NCBI_SCOPE

#endif // OBJECTS_SEQFEAT_SUBSOURCE_HPP
/* Original file checksum: lines: 94, chars: 2578, CRC32: 1c534244 */