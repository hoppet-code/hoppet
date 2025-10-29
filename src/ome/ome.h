/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

/**
 * \file
 * \brief Main header for libome
 *
 * \details
 * This header includes all necessary further headers to use the library.
 * It works for both C and C++. If it is included from C sources, only the
 * C interface will be exposed.
 */

#ifndef LIBOME_OME_H
#define LIBOME_OME_H

#ifdef __cplusplus
#include <ome/ome_type_aliases.h>
#include <ome/mellin.h>
#include <ome/integration_engine_gsl.h>
#endif

#include <ome/AqqQNSEven.h>
#include <ome/AqqQNSOdd.h>
#include <ome/AQqPS.h>
#include <ome/AqqQPS.h>
#include <ome/AqgQ.h>
#include <ome/AgqQ.h>
#include <ome/AggQ.h>
#include <ome/AQg.h>

#include <ome/polAqqQNSEven.h>
#include <ome/polAqqQNSOdd.h>
#include <ome/polAQqPS.h>
#include <ome/polAqqQPS.h>
#include <ome/polAqgQ.h>
#include <ome/polAgqQ.h>
#include <ome/polAggQ.h>
#include <ome/polAQg.h>

#endif
