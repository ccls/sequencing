/*  $Id: wsdlparser.cpp 282780 2011-05-16 16:02:27Z gouriano $
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
* Author: Andrei Gourianov
*
* File Description:
*   WSDL parser
*
* ===========================================================================
*/

#include <ncbi_pch.hpp>
#include "exceptions.hpp"
#include "wsdlparser.hpp"
#include "tokens.hpp"
#include "module.hpp"
#include "datatool.hpp"
#include <serial/error_codes.hpp>


#define NCBI_USE_ERRCODE_X   Serial_Parsers

BEGIN_NCBI_SCOPE

/////////////////////////////////////////////////////////////////////////////
// WSDL Parser

WSDLParser::WSDLParser(WSDLLexer& lexer)
    : XSDParser(lexer)
{
    m_SrcType = eWsdl;
    m_ParsingTypes = false;
    m_ParsingOutput = false;
}

WSDLParser::~WSDLParser(void)
{
}

void WSDLParser::BuildDocumentTree(CDataTypeModule& module)
{
    ParseHeader();
    CopyComments(module.Comments());

    TToken tok;
    for (;;) {
        tok = GetNextToken();
        switch ( tok ) {
        case K_TYPES:
            ParseTypes(module);
            break;
        case K_MESSAGE:
            ParseMessage();
            break;
        case K_PORTTYPE:
            CreateTypeDefinition( DTDEntity::eWsdlInterface);
            break;
        case K_BINDING:
            CreateTypeDefinition( DTDEntity::eWsdlBinding);
            break;
        case K_SERVICE:
            ParseService();
            break;
        case K_ENDOFTAG:
            m_ParsingTypes = true;
            ProcessNamedTypes();
            m_ParsingTypes = false;
            ProcessEndpointTypes();
            CollectDataObjects();
            return;
        case T_EOF:
            ParseError("Unexpected end-of-file", "keyword");
            return;
        default:
            ParseError("Invalid keyword", "keyword");
            return;
        }
    }
}

void WSDLParser::ParseHeader(void)
{
// xml header
    TToken tok = GetNextToken();
    if (tok == K_XML) {
        for ( ; tok != K_ENDOFTAG; tok=GetNextToken())
            ;
        tok = GetNextToken();
    } else {
        ERR_POST_X(4, "LINE " << Location() << " XML declaration is missing");
    }
// schema    
    if (tok != K_DEFINITIONS) {
        ParseError("Unexpected token", "definitions");
    }
    for ( tok = GetNextToken(); tok == K_ATTPAIR || tok == K_XMLNS; tok = GetNextToken()) {
        if (tok == K_ATTPAIR) {
            if (IsAttribute("targetNamespace")) {
                m_TargetNamespace = m_Value;
            }
        }
    }
    if (tok != K_CLOSING) {
        ParseError("tag closing");
    }
}

void WSDLParser::ParseTypes(CDataTypeModule& module)
{
    TToken tok = GetRawAttributeSet();
    if (tok == K_CLOSING) {
        m_ParsingTypes = true;
        while (DTDParser::GetNextToken() == K_SCHEMA) {
            WSDLLexer& l = dynamic_cast<WSDLLexer&>(Lexer());
            l.UseXSDLexer(true);
            BeginScope(NULL);
            XSDParser::BuildDocumentTree(module);
            EndScope();
            l.UseXSDLexer(false);
        }
        m_ParsingTypes = false;
        tok = GetNextToken();
        if (tok != K_ENDOFTAG) {
            ParseError("Unexpected token", "end of tag");
        }
    }
}

void WSDLParser::ParseContent(DTDElement& node)
{
    bool dounk=false;
    TToken tok;
    for ( tok=GetNextToken(); ; tok=GetNextToken()) {
        dounk=false;
        switch (tok) {
        case T_EOF:
            return;
        case K_ENDOFTAG:
            return;
        case K_PORTTYPE:
            ParsePortType(node);
            break;
        case K_BINDING:
            ParseBinding(node);
            break;
        case K_OPERATION:
            ParseOperation(node);
            break;
        case K_INPUT:
            ParseInput(node);
            break;
        case K_OUTPUT:
            ParseOutput(node);
            break;
        case K_BODY:
            ParseBody(node);
            break;
        case K_HEADER:
            ParseHeader(node);
            break;
        case K_PART:
            ParsePart(node);
            break;
        case K_PORT:
            ParsePort(node);
            break;
        case K_ADDRESS:
            ParseAddress(node);
            break;
        case K_DOCUMENTATION:
            SetCommentsIfEmpty(&(node.Comments()));
            ParseDocumentation();
            break;
        default:
            dounk = true;
            break;
        }
        if (dounk) {
            for ( tok = GetNextToken(); tok == K_ATTPAIR || tok == K_XMLNS; tok = GetNextToken())
                ;
            if (tok == K_CLOSING) {
                ParseContent(node);
            }
        }
    }
}

void WSDLParser::ParsePortType(DTDElement& node)
{
    TToken tok = GetRawAttributeSet();
    if (tok == K_CLOSING) {
        ParseContent(node);
    }
}

void WSDLParser::ParseBinding(DTDElement& node)
{
    EElementNamespace ns = GetElementNamespace(m_ElementPrefix);
    TToken tok = GetRawAttributeSet();
    if (ns == eWsdlNamespace) {
        if (!GetAttribute("type")) {
            ParseError("Binding has no type", "type");
        }
        string id = CreateEntityId(m_Value,DTDEntity::eWsdlInterface);
        if (m_MapEntity.find(id) == m_MapEntity.end()) {
            ParseError("Unresolved entity", id.c_str());
        }
        PushEntityLexer(id);
        ParseContent(node);
    } else if (ns == eSoapNamespace) {
        if (GetAttribute("style")) {
            if (m_Value != "document") {
                ParseError("Only document style binding is supported", "document");
            }
        }
    } else {
        node.SetType(DTDElement::eWsdlUnsupportedEndpoint);
    }
    if (tok == K_CLOSING) {
        ParseContent(node);
    }
}

void WSDLParser::ParseOperation(DTDElement& node)
{
    TToken tok = GetRawAttributeSet();
    if (node.GetType() == DTDElement::eWsdlEndpoint) {
        if (!GetAttribute("name")) {
            ParseError("Operation has no name", "name");
        }
        DTDElement& item = EmbeddedElement(node, m_Value,DTDElement::eWsdlOperation);
        item.SetSourceLine(Lexer().CurrentLine());
        if (tok == K_CLOSING) {
            ParseContent(item);
        }
        return;
    }
    if (GetAttribute("soapAction")) {
        string value(m_ValuePrefix);
        if (!value.empty()) {
            value += ':';
        }
        value += m_Value;
        if (!value.empty()) {
            DTDElement& action = EmbeddedElement(node, "#soapaction", DTDElement::eString);
            action.SetDefault(value);
        }
    }
    if (tok == K_CLOSING) {
        ParseContent(node);
    }
}

void WSDLParser::ParseBody(DTDElement& node)
{
    TToken tok = GetRawAttributeSet();
    if (tok == K_CLOSING) {
        SkipContent();
    }
}

void WSDLParser::ParseHeader(DTDElement& node)
{
    EElementNamespace ns = GetElementNamespace(m_ElementPrefix);
    TToken tok = GetRawAttributeSet();
    if (ns == eSoapNamespace) {
        if (GetAttribute("message")) {
            DTDElement::EType etype;
            etype = m_ParsingOutput ? DTDElement::eWsdlHeaderOutput : DTDElement::eWsdlHeaderInput;
            string item_name(CreateEmbeddedName(node, etype));
            DTDElement& item = m_MapElement[item_name];
            item.SetEmbedded();
            item.SetName(item_name);
            item.SetType(etype);
            item.SetSourceLine(Lexer().CurrentLine());
            AddElementContent(node,item_name);

            string msg_name(CreateWsdlName(m_Value,DTDElement::eWsdlMessage));
            DTDElement& msg = m_MapElement[msg_name];
            msg.SetName(msg_name);
            msg.SetType( DTDElement::eWsdlMessage);
            AddElementContent(item,msg_name);
            if (tok == K_CLOSING) {
                ParseContent(item);
            }
            return;
        }
    }
    if (tok == K_CLOSING) {
        SkipContent();
    }
}

void WSDLParser::ParseInput(DTDElement& node)
{
    TToken tok = GetRawAttributeSet();
    if (GetAttribute("message")) {
        string item_name(CreateEmbeddedName(node, DTDElement::eWsdlInput));
        DTDElement& item = m_MapElement[item_name];
        item.SetEmbedded();
        item.SetName(item_name);
        item.SetType(DTDElement::eWsdlInput);
        item.SetSourceLine(Lexer().CurrentLine());
        AddElementContent(node,item_name);

        string msg_name(CreateWsdlName(m_Value,DTDElement::eWsdlMessage));
        DTDElement& msg = m_MapElement[msg_name];
        msg.SetName(msg_name);
        msg.SetType( DTDElement::eWsdlMessage);
        AddElementContent(item,msg_name);
        if (tok == K_CLOSING) {
            ParseContent(item);
        }
        return;
    }
    if (tok == K_CLOSING) {
        ParseContent(node);
    }
}

void WSDLParser::ParseOutput(DTDElement& node)
{
    TToken tok = GetRawAttributeSet();
    if (GetAttribute("message")) {
        string item_name(CreateEmbeddedName(node, DTDElement::eWsdlOutput));
        DTDElement& item = m_MapElement[item_name];
        item.SetEmbedded();
        item.SetName(item_name);
        item.SetType(DTDElement::eWsdlOutput);
        item.SetSourceLine(Lexer().CurrentLine());
        AddElementContent(node,item_name);

        string msg_name(CreateWsdlName(m_Value,DTDElement::eWsdlMessage));
        DTDElement& msg = m_MapElement[msg_name];
        msg.SetName(msg_name);
        msg.SetType( DTDElement::eWsdlMessage);
        AddElementContent(item,msg_name);
        if (tok == K_CLOSING) {
            ParseContent(item);
        }
        return;
    }
    if (tok == K_CLOSING) {
        m_ParsingOutput = true;
        ParseContent(node);
        m_ParsingOutput = false;
    }
}

void WSDLParser::ParsePart(DTDElement& node)
{
    TToken tok = GetRawAttributeSet();
    if (GetAttribute("element")) {
        AddElementContent(node,m_Value);
    } else {
        bool ok = false;
        if (GetAttribute("name")) {
            string name(m_Value);
            if (GetAttribute("type")) {
                DTDElement& item = EmbeddedElement(node, name, DTDElement::eUnknown);
                if (DefineElementType(item)) {
                    ok = true;
                }
            }
        }
        if (!ok) {
            ParseError("Part has no element", "element");
        }
    }
    if (tok == K_CLOSING) {
        SkipContent();
    }
}

void WSDLParser::ParsePort(DTDElement& node)
{
    TToken tok = GetRawAttributeSet();
    if (!GetAttribute("name")) {
        ParseError("Port has no name", "name");
    }
    string name(CreateWsdlName(m_Value,DTDElement::eWsdlEndpoint));
    DTDElement& item = m_MapElement[name];
    item.SetName(m_Value);
    item.SetSourceLine(Lexer().CurrentLine());
    item.SetType( DTDElement::eWsdlEndpoint);
    item.SetNamespaceName(m_TargetNamespace);
    if (!GetAttribute("binding")) {
        ParseError("Port has no binding", "binding");
    }
    item.SetTypeName(CreateEntityId(m_Value,DTDEntity::eWsdlBinding));
    AddElementContent(node,name);
    if (tok == K_CLOSING) {
        ParseContent(item);
    }
}

void WSDLParser::ParseAddress(DTDElement& node)
{
    TToken tok = GetRawAttributeSet();
    if (!GetAttribute("location")) {
        ParseError("Port has no location", "location");
    }
    string value(m_ValuePrefix);
    if (!value.empty()) {
        value += ':';
    }
    value += m_Value;
    DTDElement& item = EmbeddedElement(node, "#location", DTDElement::eString);
    item.SetSourceLine(Lexer().CurrentLine());
    item.SetDefault(value);
    if (tok == K_CLOSING) {
        ParseContent(item);
    }
}

string WSDLParser::CreateWsdlName(const string& name, DTDElement::EType type)
{
    string id;
    switch (type) {
    case DTDElement::eWsdlService:
        id = string("service:") + name;
        break;
    case DTDElement::eWsdlEndpoint:
        id = string("endpoint:") + name;
        break;
    case DTDElement::eWsdlOperation:
        id = string("operation:") + name;
        break;
    case DTDElement::eWsdlHeaderInput:
        id = string("headerinput:") + name;
        break;
    case DTDElement::eWsdlInput:
        id = string("input:") + name;
        break;
    case DTDElement::eWsdlHeaderOutput:
        id = string("headeroutput:") + name;
        break;
    case DTDElement::eWsdlOutput:
        id = string("output:") + name;
        break;
    case DTDElement::eWsdlMessage:
        id = string("message:") + name;
        break;
    default:
        id = name;
        break;
    }
    return id;
}

string WSDLParser::CreateEmbeddedName(DTDElement& node, DTDElement::EType type)
{
    return CreateWsdlName( CreateTmpEmbeddedName(
        node.GetName(), node.GetContent().size()),type);
}

DTDElement& WSDLParser::EmbeddedElement(DTDElement& node,
    const string& name, DTDElement::EType type)
{
    const list<string>& cont = node.GetContent();
    ITERATE( list<string>, i, cont) {
        if (m_MapElement.find( *i) != m_MapElement.end() &&
            m_MapElement[*i].GetName() == name) {
            return m_MapElement[*i];
        }
    }
    string item_name(CreateEmbeddedName(node, type));
    DTDElement& item = m_MapElement[item_name];
    item.SetEmbedded();
    item.SetNamed();
    item.SetName(name);
    item.SetType(type);
    AddElementContent(node,item_name);
    return item;
}

void WSDLParser::ParseMessage(void)
{
    TToken tok = GetRawAttributeSet();
    if (!GetAttribute("name")) {
        ParseError("Message has no name", "name");
    }

    string msg_name(CreateWsdlName(m_Value,DTDElement::eWsdlMessage));
    DTDElement& node = m_MapElement[msg_name];
    node.SetName(msg_name);
    node.SetType( DTDElement::eWsdlMessage);
    node.SetSourceLine(Lexer().CurrentLine());
    if (tok == K_CLOSING) {
        ParseContent(node);
    }
}

void WSDLParser::ParseService(void)
{
    TToken tok = GetRawAttributeSet();
    if (!GetAttribute("name")) {
        ParseError("Service has no name", "name");
    }
//    ((CDataTool*)CNcbiApplication::Instance())->SetDefaultNamespace(string("NCBI_NS_NCBI::objects") + "::" + m_Value);
    DTDElement& node = m_MapElement[CreateWsdlName(m_Value,DTDElement::eWsdlService)];
    node.SetName(m_Value);
    node.SetSourceLine(Lexer().CurrentLine());
    node.SetType( DTDElement::eWsdlService);
    if (tok == K_CLOSING) {
        ParseContent(node);
    }
}

AbstractLexer* WSDLParser::CreateEntityLexer(
    CNcbiIstream& in, const string& name, bool autoDelete /*=true*/)
{
    WSDLEntityLexer* l = new WSDLEntityLexer(in,name);
    l->UseXSDLexer(m_ParsingTypes);
    return l;
}

void WSDLParser::ProcessEndpointTypes(void)
{
    map<string,DTDElement>::iterator i;
    bool found = false;
    for (i = m_MapElement.begin(); i != m_MapElement.end(); ++i) {
        DTDElement& node = i->second;
        if (node.GetType() == DTDElement::eWsdlEndpoint) {
            found = true;
            PushEntityLexer(node.GetTypeName());
            ParseContent(node);
        }
    }
    if (found) {
        // remove unsupported endpoints
        for (i = m_MapElement.begin(); i != m_MapElement.end(); ++i) {
            DTDElement& node = i->second;
            if (node.GetType() == DTDElement::eWsdlService) {
                list<string> refs = node.GetContent();
                ITERATE(list<string>,r,refs) {
                    DTDElement& refNode = m_MapElement[*r];
                    if (refNode.GetType() == DTDElement::eWsdlUnsupportedEndpoint) {
                        node.RemoveContent(*r);
                        m_MapElement.erase(*r);
                    }
                }
            }
        }
    } else {
        // no endpoint spec found
        // create default (?)
        map<string,DTDEntity>::iterator e;
        for (e = m_MapEntity.begin(); e != m_MapEntity.end(); ++e) {
            DTDEntity& ent = e->second;
            string ename(e->first);
            if (ent.GetType() == DTDEntity::eWsdlInterface) {
                string name(CreateWsdlName(ent.GetName(),DTDElement::eWsdlEndpoint));
                DTDElement& item = m_MapElement[name];
                item.SetName(ent.GetName());
                item.SetType( DTDElement::eWsdlEndpoint);
                item.SetNamespaceName(m_TargetNamespace);
                item.SetTypeName(ename);
                found = true;
            }
        }
        if (found) {
            ProcessEndpointTypes();
        }
    }
}

void WSDLParser::CollectDataObjects(void)
{
    map<string,DTDElement>::iterator i;
    for (i = m_MapElement.begin(); i != m_MapElement.end(); ++i) {
        DTDElement& node = i->second;
        if (node.GetType() == DTDElement::eWsdlEndpoint ||
            node.GetType() == DTDElement::eWsdlHeaderInput ||
            node.GetType() == DTDElement::eWsdlInput ||
            node.GetType() == DTDElement::eWsdlHeaderOutput ||
            node.GetType() == DTDElement::eWsdlOutput) {
            CollectDataObjects(node,node);
        }
    }
}

void WSDLParser::CollectDataObjects(DTDElement& agent, DTDElement& node)
{
    const list<string> refs = node.GetContent();
    for (list<string>::const_iterator i= refs.begin(); i != refs.end(); ++i) {
        DTDElement& refNode = m_MapElement[*i];
        if (refNode.GetType() < DTDElement::eWsdlService) {
            if (&agent != &node &&
//                !refNode.IsEmbedded() &&
//                node.GetType() != DTDElement::eWsdlOperation &&
                find(agent.GetContent().begin(), agent.GetContent().end(),*i) ==
                     agent.GetContent().end()) {
                agent.AddContent( *i);
            }
        } else {
            CollectDataObjects(agent,refNode);
        }
    }
}

END_NCBI_SCOPE