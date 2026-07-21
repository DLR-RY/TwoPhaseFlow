/*---------------------------------------------------------------------------*\
  Runtime-selection compatibility helpers for OpenFOAM v2606.

  These macros retain OpenFOAM's normal first-registration-wins behaviour,
  but do not emit a duplicate-entry stack trace when another loaded library
  has already registered the same lookup key.
\*---------------------------------------------------------------------------*/

#ifndef Foam_addIfAbsentToRunTimeSelectionTable_H
#define Foam_addIfAbsentToRunTimeSelectionTable_H

#include "addToRunTimeSelectionTable.H"


//- Add to a construction table only when the typeName key is not present
#define addIfAbsentToRunTimeSelectionTable(baseType, thisType, argNames)       \
                                                                               \
    struct addIfAbsent##thisType##argNames##ConstructorTo##baseType##Table    \
    {                                                                          \
        addIfAbsent##thisType##argNames##ConstructorTo##baseType##Table()     \
        {                                                                      \
            baseType::argNames##ConstructorTablePtr_construct(true);          \
            const ::Foam::word& key = thisType::typeName;                     \
            if (!baseType::argNames##ConstructorTablePtr_->found(key))        \
            {                                                                  \
                baseType::argNames##ConstructorTablePtr_->insert              \
                (                                                              \
                    key,                                                       \
                    baseType::add##argNames##ConstructorToTable<thisType>::New\
                );                                                             \
            }                                                                  \
        }                                                                      \
                                                                               \
        ~addIfAbsent##thisType##argNames##ConstructorTo##baseType##Table()    \
        {                                                                      \
            baseType::argNames##ConstructorTablePtr_construct(false);         \
        }                                                                      \
                                                                               \
        addIfAbsent##thisType##argNames##ConstructorTo##baseType##Table       \
        (                                                                      \
            const addIfAbsent##thisType##argNames##ConstructorTo##baseType##  \
                Table&                                                        \
        ) = delete;                                                            \
                                                                               \
        void operator=                                                         \
        (                                                                      \
            const addIfAbsent##thisType##argNames##ConstructorTo##baseType##  \
                Table&                                                        \
        ) = delete;                                                            \
    }                                                                          \
    addIfAbsent##thisType##argNames##ConstructorTo##baseType##Table_


//- Add to a construction table only when the named lookup key is not present
#define addNamedIfAbsentToRunTimeSelectionTable\
(baseType, thisType, argNames, lookupName)                                    \
                                                                               \
    struct addIfAbsent##thisType##argNames##ConstructorTo##baseType##Table_## \
        lookupName                                                            \
    {                                                                          \
        addIfAbsent##thisType##argNames##ConstructorTo##baseType##Table_##    \
            lookupName()                                                       \
        {                                                                      \
            baseType::argNames##ConstructorTablePtr_construct(true);          \
            const ::Foam::word key(#lookupName);                              \
            if (!baseType::argNames##ConstructorTablePtr_->found(key))        \
            {                                                                  \
                baseType::argNames##ConstructorTablePtr_->insert              \
                (                                                              \
                    key,                                                       \
                    baseType::add##argNames##ConstructorToTable<thisType>::New\
                );                                                             \
            }                                                                  \
        }                                                                      \
                                                                               \
        ~addIfAbsent##thisType##argNames##ConstructorTo##baseType##Table_##   \
            lookupName()                                                       \
        {                                                                      \
            baseType::argNames##ConstructorTablePtr_construct(false);         \
        }                                                                      \
                                                                               \
        addIfAbsent##thisType##argNames##ConstructorTo##baseType##Table_##    \
            lookupName                                                        \
        (                                                                      \
            const addIfAbsent##thisType##argNames##ConstructorTo##baseType##  \
                Table_##lookupName&                                            \
        ) = delete;                                                            \
                                                                               \
        void operator=                                                         \
        (                                                                      \
            const addIfAbsent##thisType##argNames##ConstructorTo##baseType##  \
                Table_##lookupName&                                            \
        ) = delete;                                                            \
    }                                                                          \
    addIfAbsent##thisType##argNames##ConstructorTo##baseType##Table_##        \
        lookupName##_


#endif

// ************************************************************************* //
