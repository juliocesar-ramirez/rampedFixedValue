/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dimensionSet.H"
#include "rampedFixedValueFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "scalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "zero.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::rampedFixedValueFvPatchField<Type>::t() const
{
    return this->db().time().timeOutputValue();
}

template<class Type>
Foam::scalar Foam::rampedFixedValueFvPatchField<Type>::coeficienteX() const
{
    scalar y= (t()-tiempoInicial_)/(tiempoFinal_-tiempoInicial_);
    return min(1,max(0,y));

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Type>
Foam::rampedFixedValueFvPatchField<Type>::rampedFixedValueFvPatchField(
    const fvPatch &p, const DimensionedField<Type, volMesh> &iF)
    : fixedValueFvPatchField<Type>(p, iF), scalarData_(0.0), tiempoFinal_(0.0),
      tiempoInicial_(0.0), data_(Zero), valorAlto_(Zero),valorBajo_(Zero), fieldData_(p.size(), Zero) {}

template<class Type>
Foam::rampedFixedValueFvPatchField<Type>::
rampedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    scalarData_(dict.lookup<scalar>("scalarData")),
    tiempoFinal_(dict.lookup<scalar>("tiempoFinal")),
    tiempoInicial_(dict.lookup<scalar>("tiempoInicial")),
    data_(dict.lookup<Type>("data")),
    valorAlto_(dict.lookup<Type>("valorAlto")),
    valorBajo_(dict.lookup<Type>("valorBajo")),
    fieldData_("fieldData", dict, p.size())
{


    fixedValueFvPatchField<Type>::evaluate();

    /*
    // Initialise with the value entry if evaluation is not possible
    fvPatchField<Type>::operator=
    (
        Field<Type>("value", dict, p.size())
    );
    */
}


template<class Type>
Foam::rampedFixedValueFvPatchField<Type>::
rampedFixedValueFvPatchField
(
    const rampedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    scalarData_(ptf.scalarData_),
    tiempoInicial_(ptf.tiempoInicial_),
    tiempoFinal_(ptf.tiempoFinal_),
    data_(ptf.data_),
    valorAlto_(ptf.valorAlto_),
    valorBajo_(ptf.valorBajo_),
    fieldData_(mapper(ptf.fieldData_))
{}

template <class Type>
Foam::rampedFixedValueFvPatchField<Type>::rampedFixedValueFvPatchField(
    const rampedFixedValueFvPatchField<Type> &ptf,
    const DimensionedField<Type, volMesh> &iF)
    : fixedValueFvPatchField<Type>(ptf, iF), scalarData_(ptf.scalarData_),
      tiempoInicial_(ptf.tiempoInicial_),tiempoFinal_(ptf.tiempoFinal_),
      data_(ptf.data_),valorAlto_(ptf.valorAlto_),valorBajo_(ptf.valorBajo_),
      fieldData_(ptf.fieldData_) {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * *
// * //

template <class Type>
void Foam::rampedFixedValueFvPatchField<Type>::autoMap(
    const fvPatchFieldMapper &m) {
  fixedValueFvPatchField<Type>::autoMap(m);
  m(fieldData_, fieldData_);
}

template<class Type>
void Foam::rampedFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const rampedFixedValueFvPatchField<Type>& tiptf =
        refCast<const rampedFixedValueFvPatchField<Type>>(ptf);

    fieldData_.rmap(tiptf.fieldData_, addr);
}


template<class Type>
void Foam::rampedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    fixedValueFvPatchField<Type>::operator==
    (
        valorBajo_
        +(valorAlto_-valorBajo_)*coeficienteX()
    );


    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::rampedFixedValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "scalarData", scalarData_);
    writeEntry(os, "tiempoInicial", tiempoInicial_);
    writeEntry(os, "tiempoFinal", tiempoFinal_);
    writeEntry(os, "data", data_);
    writeEntry(os, "valorAlto", valorAlto_);
    writeEntry(os, "valorBajo", valorBajo_);
    writeEntry(os, "fieldData", fieldData_);
    writeEntry(os, "value", *this);
}



// ************************************************************************* //
