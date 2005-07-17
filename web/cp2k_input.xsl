<?xml version="1.0" encoding="ISO-8859-1"?>

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="html" encoding="iso-8859-1" indent="yes"
 doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"
 doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>

<xsl:template match="/INPUT">
 <xsl:for-each select="SECTION">
  <html xmlns="http://www.w3.org/1999/xhtml">
   <head>
    <title>CP2K input</title>
    <style type="text/css">
     body {background-color: #eeeeee}
     ul.none {list-style-type: none}
     ul.disc {list-style-type: disc}
     ul.circle {list-style-type: circle}
    </style>
   </head>
   <body>
    <h1 align="center">CP2K input description</h1>
    <h2>Version information</h2>
     This HTML manual was generated automatically from a CP2K executable
     compiled on <xsl:value-of select="/INPUT/COMPILE_DATE"/> using the
     -xml command line option. Thus the manual describes exactly this
     version of the code. The latest CVS log file entry found was
     <xsl:value-of select="/INPUT/COMPILE_LASTCVS"/>.
    <h2>Input structure</h2>
     All sections and keywords that can be part of a CP2K input file are shown
     with their allowed nestings. A detailed description can be obtained by
     clicking on the links. The links in the detailed descriptions switch
     back to the corresponding index entry. In this way a toggling between
     the index and the detailed description is feasible.
    <h2>Index of all sections and keywords</h2>
     <xsl:call-template name="section_index"></xsl:call-template>
     Back to the <a href="http://cp2k.berlios.de">CP2K home page</a>
    <hr/>
    <h2>Detailed description of all sections and keywords</h2>
     <xsl:call-template name="describe_sections"></xsl:call-template>
    <hr/>
    Back to the <a href="http://cp2k.berlios.de">CP2K home page</a>
   </body>
  </html>
 </xsl:for-each>
</xsl:template>

<xsl:template name="section_index">
 <xsl:choose>
  <xsl:when test="count(SECTION) = 0">
   <br/>
  </xsl:when>
  <xsl:otherwise>
   <xsl:for-each select="SECTION">
    <xsl:sort select="NAME"/>
    <ul class="disc">
     <li><a href="#sec_des_{generate-id(NAME)}" id="sec_ind_{generate-id(NAME)}">&amp;<xsl:value-of select="NAME"/></a>
      <br/>
      <xsl:call-template name="keyword_index"></xsl:call-template>
      <xsl:call-template name="section_index"></xsl:call-template>
     </li>
    </ul>
   </xsl:for-each>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

<xsl:template name="keyword_index">
 <xsl:for-each select="KEYWORD">
  <xsl:sort select="NAME[@type='default']"/>
  <ul class="circle">
   <li><a href="#key_des_{generate-id(NAME[@type='default'])}" id="key_ind_{generate-id(NAME[@type='default'])}"><xsl:value-of select="NAME[@type='default']"/></a></li>
  </ul>
 </xsl:for-each>
 <br/>
</xsl:template>

<xsl:template name="describe_sections">
 <xsl:choose>
  <xsl:when test="count(SECTION) = 0">
  </xsl:when>
  <xsl:otherwise>
   <xsl:for-each select="SECTION">
    <xsl:sort select="NAME"/>
    <br/>
    <a href="#sec_ind_{generate-id(NAME)}" id="sec_des_{generate-id(NAME)}"><h3>Section &amp;<xsl:value-of select="NAME"/></h3></a>
    <ul class="none">
     <li>
      <em><xsl:value-of select="DESCRIPTION"/></em>
      <xsl:call-template name="describe_keywords"></xsl:call-template>
     </li>
    </ul>
    <xsl:call-template name="describe_sections"></xsl:call-template>
   </xsl:for-each>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

<xsl:template name="describe_keywords">
 <xsl:choose>
  <xsl:when test="count(KEYWORD) = 0">
  </xsl:when>
  <xsl:otherwise>
   <h4>Keywords:</h4>
   <ul class="none">
   <li>
   <xsl:for-each select="KEYWORD">
    <xsl:sort select="NAME"/>
    <dl>
     <xsl:if test="position() != 1"><br/></xsl:if>
     <dt>
      <a href="#key_ind_{generate-id(NAME[@type='default'])}" id="key_des_{generate-id(NAME[@type='default'])}"><xsl:value-of select="NAME[@type='default']"/></a>
      <xsl:if test="NAME[@type='alias']">
       (alias: <xsl:value-of select="NAME[@type='alias']"/>)
      </xsl:if>
     </dt>
     <dd><p><em><xsl:value-of select="DESCRIPTION"/></em></p></dd>
     <dd>
      <p>Data type: <big><tt><xsl:value-of select="DATA_TYPE/@kind"/></tt></big>
      <xsl:if test="DATA_TYPE/ENUMERATION">
       <p>List of valid keys:</p>
       <ul class="none">
       <xsl:for-each select="DATA_TYPE/ENUMERATION/ITEM">
        <xsl:sort select="NAME"/>
        <dl>
         <dt><big><tt><xsl:value-of select="NAME"/></tt></big></dt>
         <dd><em><xsl:value-of select="DESCRIPTION"/></em></dd>
        </dl>
       </xsl:for-each>
       </ul>
      </xsl:if>
      </p>
     </dd>
     <dd><p>Usage example: <big><tt><xsl:value-of select="USAGE"/></tt></big></p></dd>
     <xsl:if test="string-length(DEFAULT_VALUE) > 0">
      <dd><p>Default value: <big><tt><xsl:value-of select="DEFAULT_VALUE"/></tt></big></p></dd>
     </xsl:if>
    </dl>
   </xsl:for-each>
   </li>
   </ul>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

</xsl:stylesheet>
