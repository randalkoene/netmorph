/*
  © Copyright 2024 Randal A. Koene <randalk@netmorph.org>
  
  With design assistance from J. van Pelt & A. van Ooyen, and support
  from the Netherlands Organization for Scientific Research (NWO)
  Program Computational Life Sciences grant CLS2003 (635.100.005) and
  from the EC Marie Curie Research and Training Network (RTN)
  NEURoVERS-it 019247.

  This file is part of NETMORPH.

  NETMORPH is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  NETMORPH is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with NETMORPH.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include <memory>

/**
 * Glue-class to seamlessly pass logging messages from Netmorph
 * to the standard NES Logger.
 */
class Netmorph2NESLogging {
public:
  virtual void error(const std::string & msg) = 0;
  virtual void warning(const std::string & msg) = 0;
  virtual void report(const std::string & msg) = 0;
  virtual void progress(const std::string & msg) = 0;
};

extern std::unique_ptr<Netmorph2NESLogging> embedlog; /** Used if Netmorph should behave as an embedded component. */
