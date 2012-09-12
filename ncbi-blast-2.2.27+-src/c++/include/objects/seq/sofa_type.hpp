/*  $Id: sofa_type.hpp 371868 2012-08-13 15:10:25Z rafanovi $
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
 * Authors:  Frank Ludwig
 *
 * File Description:  Sequence Ontology Type Mapping
 *
 */

#ifndef OBJECTS___SOFATYPE__HPP
#define OBJECTS___SOFATYPE__HPP

#include <objects/seqfeat/SeqFeatData.hpp>

BEGIN_NCBI_SCOPE
BEGIN_objects_SCOPE

//  ============================================================================
struct NCBI_SEQ_EXPORT SofaType
//  ============================================================================
{
    unsigned int m_id;
    string m_name;

    SofaType() :
    m_id( 0 ),
    m_name( "invalid" )
    {};

    SofaType(
        unsigned int id,
        const string& name ) :
    m_id( id ),
    m_name( name )
    {};
};

END_objects_SCOPE
END_NCBI_SCOPE

#endif  // OBJECTS___SOFATYPE__HPP