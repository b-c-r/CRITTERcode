################################################################################
#   runs mod09h_rrpe_fit                                                       #
#                                                                              #
#    Copyright (C) 2025 Björn C. Rall (https://orcid.org/0000-0002-3191-8389)  #
#                                                                              #
#    This program is free software: you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License as published by      #
#    the Free Software Foundation, either version 3 of the License, or         #
#    (at your option) any later version.                                       #
#                                                                              #
#    This program is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#    GNU General Public License for more details.                              #
#                                                                              #
#    You should have received a copy of the GNU General Public License         #
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.    #
################################################################################

mod09h_fit_ie <- mod09h_rrpe_fit(
  n_eaten = fr_data_ie$n_eaten,                                                  # data: number of prey eaten, as integer
  n_initial = fr_data_ie$n_initial,                                              # data: number of prey provided initially, as integer
  n_rings = fr_data_ie$ring_count,                                               # data: number of habitat rings provided as structure
  complexity = fr_data_ie$complexity_level                                       # data: complexity levels
)
