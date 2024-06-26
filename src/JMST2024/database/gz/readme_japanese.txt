V/50           Bright Star Catalogue, 5th Revised Ed.     (Hoffleit+, 1991)
================================================================================
The Bright Star Catalogue,  5th Revised Ed. (Preliminary Version)
     Hoffleit D., Warren Jr W.H.
    <Astronomical Data Center, NSSDC/ADC (1991)>
    =1964BS....C......0H
================================================================================
ADC_Keywords: Combined data ; Stars, bright


Description (prepared by Wayne H. Warren Jr., 1991 June 28):

    明るい星カタログ(BSC)は、6.5等星より明るい星の天文・天体物理学の
    基礎データ源として広く利用されています。このカタログには、広く
    使われている他のカタログに含まれている星の同定、二重星・多重星の同定、
    変光星の表示と変光星の同定、B1900.0 と J2000.0 の赤道上の位置が
    含まれています。 6、銀河座標、UBVRI光電測光データ（ある場合）、
    Morgan-Keenan (MK) 分類によるスペクトル型、固有運動 $(\mathrm{J} 2000.0)$ 、
    視差、半径・回転速度データ、既知の非単一星に対する多重星情報（成分数、
    分離、等級差）などが含まれています。データファイルの他に、個々のエントリ
    についてのより詳細な情報を提供する備考ファイルも充実しています。この情報
    には、星の名前、色、スペクトル、変動の詳細、連星の特徴、伴星の半径と回転速度、
    重複情報、動的視差、星の寸法（半径と直径）、偏光、星群や星団への所属などが
    含まれています。備考の有無はメインデータファイルでフラグが立てられる。

    BSCは9110天体を含み、そのうち9096天体は恒星である（1908年のオリジナル編集で
    カタログ化された14天体は新星や銀河系外天体で、番号付けを保存するために
    残されているが、そのデータはほとんど省略されている）。また、備考欄は
    メインカタログよりわずかに多い。今回の編纂では、多くの新データが追加され、
    備考欄もかなり拡大された。

    この明るい星カタログ第5版の暫定版は、Hoffleit (1982, Yale University
    Observatory) の出版版および機械読み取り可能な版に取って代わり、この版の最終版が
    完成するまでの使用を意図したものである。天文データセンターのCD-ROMで配布するためにのみ利用可能です。

    このカタログの簡単な書式説明は、暫定版のみに適用されます。最終版ではフォーマットが変更される予定です。


Author's addresses:
    Dorrit Hoffleit
        Department of Astronomy
        Yale University
    Wayne H. Warren Jr.
        ST Systems Corporation
        National Space Science Data Center
        NASA Goddard Space Flight Center


File Summary:
--------------------------------------------------------------------------------
 FileName    Lrecl    Records    Explanations
--------------------------------------------------------------------------------
ReadMe          80          .    This file
catalog        197       9110    The main part of the Catalogue
notes          132       9190    Remarks
--------------------------------------------------------------------------------

See also:
    V/36 : Supplement to the Bright Star Catalogue  (Hoffleit+ 1983)

Byte-by-byte Description of file: catalog
--------------------------------------------------------------------------------
   Bytes Format  Units   Label    Explanations
--------------------------------------------------------------------------------
   1-  4  I4     ---     HR       番号, Harvard Revised Number, Bright Star Number
   5- 14  A10    ---     Name     名前, generally Bayer and/or Flamsteed name
  15- 25  A11    ---     DM       カタログの一つのDurchmusterungの識別, Durchmusterung Identification (zone in bytes 17-19), 
  26- 31  I6     ---     HD       別カタログの番号, [1/225300]? Henry Draper Catalog Number
  32- 37  I6     ---     SAO      別カタログの番号, [1/258997]? SAO Catalog Number
  38- 41  I4     ---     FK5      別カタログの番号, ? FK5 star Number
      42  A1     ---     IRflag   赤外光源かどうか, [I] I if infrared source
      43  A1     ---   r_IRflag  *[ ':] Coded reference for infrared source
      44  A1     ---    Multiple *多重星かどうか, [AWDIRS] Double or multiple-star code
  45- 49  A5     ---     ADS      二重星カタログの指定, Aitken's Double Star Catalog (ADS) designation
  50- 51  A2     ---     ADScomp  ADS number components
  52- 60  A9     ---     VarID    変光星, Variable star identification
  61- 62  I2     h       RAh1900  ?Hours RA, equinox B1900, epoch 1900.0 (1)
  63- 64  I2     min     RAm1900  ?Minutes RA, equinox B1900, epoch 1900.0 (1)
  65- 68  F4.1   s       RAs1900  ?Seconds RA, equinox B1900, epoch 1900.0 (1)
      69  A1     ---     DE-1900  ?Sign Dec, equinox B1900, epoch 1900.0 (1)
  70- 71  I2     deg     DEd1900  ?Degrees Dec, equinox B1900, epoch 1900.0 (1)
  72- 73  I2     arcmin  DEm1900  ?Minutes Dec, equinox B1900, epoch 1900.0 (1)
  74- 75  I2     arcsec  DEs1900  ?Seconds Dec, equinox B1900, epoch 1900.0 (1)
  76- 77  I2     h       RAh      ?Hours RA, equinox J2000, epoch 2000.0 (1)
  78- 79  I2     min     RAm      ?Minutes RA, equinox J2000, epoch 2000.0 (1)
  80- 83  F4.1   s       RAs      ?Seconds RA, equinox J2000, epoch 2000.0 (1)
      84  A1     ---     DE-      ?Sign Dec, equinox J2000, epoch 2000.0 (1)
  85- 86  I2     deg     DEd      ?Degrees Dec, equinox J2000, epoch 2000.0 (1)
  87- 88  I2     arcmin  DEm      ?Minutes Dec, equinox J2000, epoch 2000.0 (1)
  89- 90  I2     arcsec  DEs      ?Seconds Dec, equinox J2000, epoch 2000.0 (1)
  91- 96  F6.2   deg     GLON     ?Galactic longitude (1)
  97-102  F6.2   deg     GLAT     ?Galactic latitude (1)
 103-107  F5.2   mag     Vmag     ?Visual magnitude (1)
     108  A1     ---   n_Vmag    *[ HR] Visual magnitude code
     109  A1     ---   u_Vmag     [ :?] Uncertainty flag on V
 110-114  F5.2   mag     B-V      ? B-V color in the UBV system
     115  A1     ---   u_B-V      [ :?] Uncertainty flag on B-V
 116-120  F5.2   mag     U-B      ? U-B color in the UBV system
     121  A1     ---   u_U-B      [ :?] Uncertainty flag on U-B
 122-126  F5.2   mag     R-I      ? R-I   in system specified by n_R-I
     127  A1     ---   n_R-I      [CE:?D] Code for R-I system (Cousin, Eggen)
 128-147  A20    ---     SpType   Spectral type
     148  A1     ---   n_SpType   [evt] Spectral type code
 149-154  F6.3 arcsec/yr pmRA    *?Annual proper motion in RA J2000, FK5 system
 155-160  F6.3 arcsec/yr pmDE     ?Annual proper motion in Dec J2000, FK5 system
     161  A1     ---   n_Parallax [D] D indicates a dynamical parallax,
                                    otherwise a trigonometric parallax
 162-166  F5.3   arcsec  Parallax ? Trigonometric parallax (unless n_Parallax)
 167-170  I4     km/s    RadVel   ? Heliocentric Radial Velocity
 171-174  A4     ---   n_RadVel  *[V?SB123O ] Radial velocity comments
 175-176  A2     ---   l_RotVel   [<=> ] Rotational velocity limit characters
 177-179  I3     km/s    RotVel   ? Rotational velocity, v sin i
     180  A1     ---   u_RotVel   [ :v] uncertainty and variability flag on
                                    RotVel
 181-184  F4.1   mag     Dmag     ? Magnitude difference of double,
                                    or brightest multiple
 185-190  F6.1   arcsec  Sep      ? Separation of components in Dmag
                                    if occultation binary.
 191-194  A4     ---     MultID   Identifications of components in Dmag
 195-196  I2     ---     MultCnt  ? Number of components assigned to a multiple
     197  A1     ---     NoteFlag [*] a star indicates that there is a note
                                    (see file notes)
--------------------------------------------------------------------------------
Note (1): These fields are all blanks for stars removed from
    the Bright Star Catalogue (see notes).
Note on r_IRflag:
  Blank if from NASA merged Infrared Catalogue, Schmitz et al., 1978;
      ' if from Engles et al. 1982
      : if uncertain identification
Note on Multiple:
      A = Astrometric binary
      D = Duplicity discovered by occultation;
      I = Innes, Southern Double Star Catalogue (1927)
      R = Rossiter, Michigan Publ. 9, 1955
      S = Duplicity discovered by speckle interferometry.
      W = Worley (1978) update of the IDS;
Note on n_Vmag:
  blank = V on UBV Johnson system;
      R = HR magnitudes reduced to the UBV system;
      H = original HR magnitude.
Note on pmRA:
     As usually assumed, the proper motion in RA is the projected
     motion (cos(DE).d(RA)/dt), i.e. the total proper motion is
     sqrt(pmRA^2^+pmDE^2^)
Note on n_RadVel:
     V  = variable radial velocity;
     V? = suspected variable radial velocity;
     SB, SB1, SB2, SB3 = spectroscopic binaries,
                         single, double or triple lined spectra;
      O = orbital data available.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: notes
--------------------------------------------------------------------------------
   Bytes Format  Units  Label     Explanations
--------------------------------------------------------------------------------
   2-  5  I4     ---    HR        [1/9110]+= Harvard Revised (HR)
   6-  7  I2     ---    Count     Note counter (sequential for a star)
   8- 11  A4     ---    Category *[A-Z: ] Remark category abbreviation:
  13-132  A120   ---    Remark    Remarks in free form text
--------------------------------------------------------------------------------
Note on Category: the following abbreviations are used:
    C   - Colors;
    D   - Double and multiple stars;
    DYN - Dynamical parallaxes;
    G   - Group membership;
    M   - Miscellaneous.
    N   - Star names;
    P   - Polarization;
    R   - Stellar radii or diameters;
    RV  - Radial and/or rotational velocities;
    S   - Spectra;
    SB  - Spectroscopic binaries;
    VAR - Variability;
    The category abbreviation is always followed by a colon (:).
--------------------------------------------------------------------------------

Historical Notes:
  * 02-Oct-1993 at CDS (Francois Ochsenbein)
    A few corrections have been inserted from the CD-ROM version
    "Selected Astronomical Catalogs, Volume 1, 1991, directory
    /combined/bsc5 at CDS with the agreement of Wayne H. Warren Jr:
    1. The spectral type for HR 6397 is from Walborn and contained octal 032
       (control-Z) characters instead of square brackets around the "n".
    2. Two remarks have been added for 6985 and 8817
    3. Byte 197 (NoteFlag) of "catalog" file corrected for stars
         202 7126 7482 7614 8982 (removed asterisk)
         285  342  841  843  991 1181 1553 1652 2269 2271 (added asterisk)
        2837 3133 3962 4522 4789 6692 7076 7328 8306 8667 (added asterisk)
  * 02-Nov-1995 at CDS (Francois Ochsenbein):
    Documentation slightly changed to accommodate to standards, and
    two lines which were inverted in "notes" have been replaced.
================================================================================
(End)                                Francois Ochsenbein     [CDS]   02-Nov-1995
