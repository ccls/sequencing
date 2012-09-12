/* $Id: SeqTable_column_info.hpp 371868 2012-08-13 15:10:25Z rafanovi $
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

/// @file SeqTable_column_info.hpp
/// User-defined methods of the data storage class.
///
/// This file was originally generated by application DATATOOL
/// using the following specifications:
/// 'seqtable.asn'.
///
/// New methods or data members can be added to it if needed.
/// See also: SeqTable_column_info_.hpp


#ifndef OBJECTS_SEQTABLE_SEQTABLE_COLUMN_INFO_HPP
#define OBJECTS_SEQTABLE_SEQTABLE_COLUMN_INFO_HPP


// generated includes
#include <objects/seqtable/SeqTable_column_info_.hpp>

// generated classes

BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::

/////////////////////////////////////////////////////////////////////////////
class NCBI_SEQ_EXPORT CSeqTable_column_info : public CSeqTable_column_info_Base
{
    typedef CSeqTable_column_info_Base Tparent;
public:
    // constructor
    CSeqTable_column_info(void);
    // destructor
    ~CSeqTable_column_info(void);

    // Maps field id to name, returns empty string for unknown ids
    static const char* GetNameForId(int id);
    // Maps field name to id, returns -1 for unknown names
    static int GetIdForName(const string& name);

private:
    // Prohibit copy constructor and assignment operator
    CSeqTable_column_info(const CSeqTable_column_info& value);
    CSeqTable_column_info& operator=(const CSeqTable_column_info& value);

};

/////////////////// CSeqTable_column_info inline methods

// constructor
inline
CSeqTable_column_info::CSeqTable_column_info(void)
{
}


/////////////////// end of CSeqTable_column_info inline methods


END_objects_SCOPE // namespace ncbi::objects::

END_NCBI_SCOPE


#endif // OBJECTS_SEQTABLE_SEQTABLE_COLUMN_INFO_HPP