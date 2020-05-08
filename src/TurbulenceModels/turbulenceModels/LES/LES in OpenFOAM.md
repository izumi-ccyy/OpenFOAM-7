# Large Eddy Simulation in OpenFOAM

## LESeddyViscosity

![inheritancec diagram for LESeddyViscosity](https://cpp.openfoam.org/v7/classFoam_1_1LESModels_1_1LESeddyViscosity__inherit__graph.png)

### LESeddyViscosity.H

```cpp
#include "LESModel.H"
#include "eddyViscosity.H"
```

define `Ce_` and `epsilon`

```cpp
dimensionedScalar Ce_;
```

```cpp
virtual tmp<volScalarField> epsilon() const;
```

### LESeddyViscosity.C

```cpp
Ce_
(
    dimensioned<scalar>::lookupOrAddToDict
    (
        "Ce",
        this->coeffDict_,
        1.048
    )
)
```

$$
C_e = 1.048
$$

```cpp
template<class BasicTurbulenceModel>
tmp<volScalarField> LESeddyViscosity<BasicTurbulenceModel>::epsilon() const
{
    tmp<volScalarField> tk(this->k());

    return volScalarField::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        Ce_*tk()*sqrt(tk())/this->delta()
    );
}
```

$$
\epsilon = \frac{C_e k \sqrt{k}}{\Delta}
$$

## kEqn

```cpp
One equation eddy-viscosity model

Eddy viscosity SGS model using a modeled balance equation to simulate the
behaviour of k.

Reference:
\verbatim
    Yoshizawa, A. (1986).
    Statistical theory for compressible turbulent shear flows,
    with the application to subgrid modeling.
    Physics of Fluids (1958-1988), 29(7), 2152-2164.
\endverbatim

The default model coefficients are
\verbatim
    kEqnCoeffs
    {
        Ck                  0.094;
        Ce                  1.048;
    }
\endverbatim
```

### kEqn.H

```cpp
#include "LESeddyViscosity.H"
```

define `Ck_`

```cpp
dimensionedScalar Ck_;
```

define function `k`, `epsilon` and `DkEff`

```cpp
//- Return SGS kinetic energy
virtual tmp<volScalarField> k() const
{
    return k_;
}

//- Return sub-grid disipation rate
virtual tmp<volScalarField> epsilon() const;

//- Return the effective diffusivity for k
tmp<volScalarField> DkEff() const
{
    return volScalarField::New
    (
        "DkEff",
        this->nut_ + this->nu()
    );
}
```

$$
D_{k, Eff} = \nu + \nu_t
$$

### kEqn.C

```cpp
template<class BasicTurbulenceModel>
void kEqn<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Ck_*sqrt(k_)*this->delta();
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}
```

$$
\nu_t = \frac{C_k \sqrt{k}}{\Delta}
$$

```cpp
template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEqn<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}
```

$$
S_k = k
$$

```cpp
k_
(
    IOobject
    (
        IOobject::groupName("k", this->alphaRhoPhi_.group()),
        this->runTime_.timeName(),
        this->mesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    this->mesh_
),

Ck_
(
    dimensioned<scalar>::lookupOrAddToDict
    (
        "Ck",
        this->coeffDict_,
        0.094
    )
)
```

$$
C_k = 0.094
$$

```cpp
template<class BasicTurbulenceModel>
tmp<volScalarField> kEqn<BasicTurbulenceModel>::epsilon() const
{
    return volScalarField::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        this->Ce_*k()*sqrt(k())/this->delta()
    );
}
```

$$
\epsilon = \frac{C_e k \sqrt{k}}{\Delta}
$$

```cpp
template<class BasicTurbulenceModel>
void kEqn<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU(fvc::grad(U));
    volScalarField G(this->GName(), nut*(tgradU() && dev(twoSymm(tgradU()))));
    tgradU.clear();

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(this->Ce_*alpha*rho*sqrt(k_)/this->delta(), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}
```

$$
divU = \nabla \cdot \mathbf{U}
$$

$$
tgradU = \nabla \mathbf{U}
$$

$$
dev(\mathbf T) = \mathbf T - \frac{1}{3}tr(\mathbf T) \mathbf I
$$

$$
twoSymm(\mathbf T) = \mathbf T + \mathbf T^T
$$

$$
\begin{aligned}
G =& \nu_t \left(\nabla \mathbf{U} : \left(
\left(\nabla \mathbf{U} + \nabla \mathbf{U^T}\right)
-\frac{1}{3}tr
\left(\nabla \mathbf{U} + \nabla \mathbf{U^T}\right) \mathbf{I}
\right)\right)  \\
=& \nu_t \left( \nabla \mathbf{U} : \left(
\left(\nabla \mathbf{U} + \nabla \mathbf{U^T}\right)
-\frac{2}{3}
\left(\nabla \cdot \mathbf{U} \right) \mathbf{I}
\right) \right)  \\
= & \nu_t \left(\nabla \mathbf{U} : \left(\nabla \mathbf{U} + \nabla \mathbf{U^T}\right) \right)
-\frac{2}{3} \nu_t
\left(\nabla \mathbf{U} : \left(\nabla \cdot \mathbf{U} \right) \mathbf{I}\right)
\end{aligned}
$$

$$
\frac{\partial (\rho k)}{\partial t} + \nabla \cdot (\rho \mathbf{U} k) - \nabla \cdot (\rho D_{k, Eff} \nabla k) = \rho G - \frac{2}{3} \rho \nabla \cdot \mathbf{U} k - \frac{C_e \rho \sqrt{k}}{\Delta} k + k
$$

## The Smagorinsky SGS model

Description

    The Smagorinsky SGS model.

    Reference:
    \verbatim
        Smagorinsky, J. (1963).
        General circulation experiments with the primitive equations: I.
        The basic experiment*.
        Monthly weather review, 91(3), 99-164.
    \endverbatim

    The form of the Smagorinsky model implemented is obtained from the
    k-equation model assuming local equilibrium which provides estimates of both
    k and epsilon separate from the sub-grid scale viscosity:

    \verbatim
        B = 2/3*k*I - 2*nuSgs*dev(D)

    where

        D = symm(grad(U));
        k from D:B + Ce*k^3/2/delta = 0
        nuSgs = Ck*sqrt(k)*delta
    \endverbatim

    The default model coefficients are
    \verbatim
        SmagorinskyCoeffs
        {
            Ck                  0.094;
            Ce                  1.048;
        }
    \endverbatim

### Smagorinsky.H

```cpp
#include "LESModel.H"
#include "LESeddyViscosity.H"
```

define `Ck_`

```cpp
dimensionedScalar Ck_;
```

Return SGS kinetic energy

```cpp
//- Return SGS kinetic energy
virtual tmp<volScalarField> k() const
{
    return k(fvc::grad(this->U_));
}
```

### Smagorinsky.C

```cpp
#include "Smagorinsky.H"
#include "fvOptions.H"
```

defination of `k`

```cpp
template<class BasicTurbulenceModel>
tmp<volScalarField> Smagorinsky<BasicTurbulenceModel>::k
(
    const tmp<volTensorField>& gradU
) const
{
    volSymmTensorField D(symm(gradU));

    volScalarField a(this->Ce_/this->delta());
    volScalarField b((2.0/3.0)*tr(D));
    volScalarField c(2*Ck_*this->delta()*(dev(D) && D));

    return volScalarField::New
    (
        IOobject::groupName("k", this->alphaRhoPhi_.group()),
        sqr((-b + sqrt(sqr(b) + 4*a*c))/(2*a))
    );
}
```

$$
symm(\mathbf{T}) = \frac{1}{2} (T + T^T)
$$

$$
\mathbf{D} = \frac{1}{2} (\nabla \mathbf {U} + \nabla \mathbf {U}^T)
$$

$$
a = \frac{C_e}{\Delta}
$$

$$
b = \frac{2}{3} \nabla \cdot \mathbf{U}
$$

$$
dev(\mathbf D) = \mathbf D - \frac{1}{3}tr(\mathbf D) \mathbf I = \frac{1}{2} (\nabla \mathbf {U} + \nabla \mathbf {U}^T) - \frac{1}{3} \nabla \cdot \mathbf{U} \mathbf I
$$

$$
c = 2 C_k \Delta \left(dev(\mathbf{D}) : \mathbf{D} \right)
$$

$$
k = \left(\frac{-b + \sqrt{b^2 + 4ac}}{2a}\right)^2
$$

Correct $\nu_t$

```cpp
template<class BasicTurbulenceModel>
void Smagorinsky<BasicTurbulenceModel>::correctNut()
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    this->nut_ = Ck_*this->delta()*sqrt(k);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}
```

$$
\nu_t = C_k \Delta \sqrt{k}
$$

Calculate $\epsilon$

```cpp
template<class BasicTurbulenceModel>
tmp<volScalarField> Smagorinsky<BasicTurbulenceModel>::epsilon() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    return volScalarField::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        this->Ce_*k*sqrt(k)/this->delta()
    );
}
```

$$
\epsilon = \frac{C_e k \sqrt{k}}{\Delta}
$$
