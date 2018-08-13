//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_ParticleSimulationManager.cpp
//! \author Alex Robinson
//! \brief  Particle simulation manager class definition
//!
//---------------------------------------------------------------------------//

// FRENSIE Includes
#include "MonteCarlo_ParticleSimulationManager.hpp"
#include "MonteCarlo_ParticleSimulationManagerFactory.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_OpenMPProperties.hpp"
#include "Utility_LoggingMacros.hpp"
#include "Utility_DesignByContract.hpp"

namespace MonteCarlo{

// Constructor
ParticleSimulationManager::ParticleSimulationManager(
                 const std::string& simulation_name,
                 const std::string& archive_type,
                 const std::shared_ptr<const FilledGeometryModel>& model,
                 const std::shared_ptr<ParticleSource>& source,
                 const std::shared_ptr<EventHandler>& event_handler,
                 const std::shared_ptr<const WeightWindows> weight_windows,
                 const std::shared_ptr<const CollisionForcer> collision_forcer,
                 const std::shared_ptr<const SimulationProperties>& properties,
                 const uint64_t next_history,
                 const uint64_t rendezvous_number )
  : d_simulation_name( simulation_name ),
    d_archive_type( archive_type ),
    d_model( model ),
    d_collision_kernel( new CollisionKernel( model, *properties ) ),
    d_transport_kernel( new TransportKernel( model ) ),
    d_source( source ),
    d_event_handler( event_handler ),
    d_weight_windows( weight_windows ),
    d_collision_forcer( collision_forcer ),
    d_properties( properties ),
    d_next_history( next_history ),
    d_rendezvous_number( rendezvous_number ),
    d_rendezvous_batch_size( 0 ),
    d_batch_size( 0 ),
    d_end_simulation( false );
{
  // Make sure that the simulation name is valid
  testPrecondition( simulation_name.size() > 0 );
  // Make sure that the archive type is valid
  testPrecondition( archive_type.size() > 0 );
  // Make sure that the model pointer is valid
  testPrecondition( model.get() );
  // Make sure that the source pointer is valid
  testPrecondition( source.get() );
  // Make sure that the event handler is valid
  testPrecondition( event_handler.get() );
  // Make sure that the weight windows pointer is valid
  testPrecondition( weight_windows.get() );
  // Make sure that the collision forcer pointer is valid
  testPrecondition( collision_forcer.get() );
  // Make sure that the properties pointer is valid
  testPrecondition( properties.get() );

  // Calculate the rendezvous batch size
  uint64_t number_of_histories = d_properties->getNumberOfHistories();

  if( number_of_histories > 0 )
  {
    d_rendezvous_batch_size =
      number_of_histories/d_properties->getMinNumberOfRendezvous();

    if( d_rendezvous_batch_size > d_properties->getMaxRendezvousBatchSize() )
      d_rendezvous_batch_size = d_properties->getMaxRendezvousBatchSize();
  }
  else
    d_rendezvous_batch_size = d_properties->getMaxRendezvousBatchSize();

  // Calculate the batch size
  d_batch_size =
    d_rendezvous_batch_size/d_properties->getMinNumberOfBatchesPerRendezvous();

  if( d_batch_size > d_properties->getMaxBatchSize() )
    d_batch_size = d_properties->getMaxBatchSize();
}

// Return the rendezvous batch size
uint64_t ParticleSimulationManager::getRendezvousBatchSize() const
{
  return d_rendezvous_batch_size;
}

// Return the batch size
uint64_t ParticleSimulationManager::getBatchSize() const
{
  return d_batch_size;
}

// Set the batch size
void ParticleSimulationManager::setBatchSize( const uint64_t batch_size )
{
  d_batch_size = batch_size;
}

// Increment the next history
void ParticleSimulationManager::incrementNextHistory( const uint64_t increment_size )
{
  d_next_history += increment_size;
}
  
// Return the model
const FilledGeometryModel& ParticleSimulationManager::getModel() const
{
  return *d_model;
}

// Return the source
const ParticleSource& ParticleSimulationManager::getSource() const
{
  return *d_source;
}

// Return the event handler
const EventHandler& ParticleSimulationManager::getEventHandler() const
{
  return *d_event_handler;
}

// Return the event handler
EventHandler& ParticleSimulationManager::getEventHandler()
{
  return *d_event_handler;
}

// Return the next history that will be completed
uint64_t ParticleSimulationManager::getNextHistory() const
{
  return d_next_history;
}

// Return the number of rendezvous
uint64_t ParticleSimulationManager::getNumberOfRendezvous() const
{
  return d_rendezvous_number;
}

// Return the rendezvous batch size
uint64_t ParticleSimulationManager::getRendezvousBatchSize() const
{
  return d_rendezvous_batch_size;
}

// Return the batch size
uint64_t ParticleSimulationManager::getBatchSize() const
{
  return d_batch_size;
}

// Run the simulation set up by the user
void ParticleSimulationManager::runSimulation()
{
  FRENSIE_LOG_NOTIFICATION( "Simulation started. " );
  FRENSIE_FLUSH_ALL_LOGS();

  // Set up the random number generator for the number of threads requested
  Utility::RandomNumberGenerator::createStreams();

  // Enable source thread support
  d_source->enableThreadSupport( Utility::OpenMPProperties::getRequestedNumberOfThreads() );

  // Enable event handler thread support
  d_event_handler->enableThreadSupport( Utility::OpenMPProperties::getRequestedNumberOfThreads() );

  // Conduct the first rendezvous (for caching only)
  this->rendezvous();

  // The simulation has started
  d_event_handler->updateObserversFromParticleSimulationStartedEvent();

  uint64_t next_rendezvous_history = d_next_history + d_rendezvous_batch_size;
  
  while( !d_event_handler->isSimulationComplete() )
  {
    if( !d_end_simulation )
    {
      this->runSimulationBatch( d_next_history, d_next_history+d_batch_size );

      this->incrementNextHistory( d_batch_size );
    }
    else // end the simulation if requested (from signal handler)
      break;
    
    if( next_rendezvous_history <= d_next_history )
    {
      this->rendezvous();

      next_rendezvous_history += d_rendezvous_batch_size;
    }
  }

  // The simulation has finished
  d_event_handler->updateObserversFromParticleSimiulationStoppedEvent();

  if( !d_end_simulation )
  {
    FRENSIE_LOG_NOTIFICATION( "Simulation finished. " );
  }
  else
  {
    FRENSIE_LOG_WARNING( "Simulation terminated. " );
  }
  FRENSIE_FLUSH_ALL_LOGS();

  this->rendezvous();
}

// Reduce distributed data
void ParticleSimulationManager::reduceData( const Utility::communicator& comm,
                                            const int root_process )
{
  comm.barrier();
  
  d_source->reduceData( comm, root_process );
  d_event_handler->reduceObserverData( comm, root_process );

  comm.barrier();
}

// Rendezvous (cache state)
void ParticleSimulationManager::rendezvous()
{
  std::string archive_name( d_simulation_name );
  archive_name += "_rendezvous_";
  archive_name += Utility::toString( d_rendezvous_number );
  archive_name += ".";
  archive_name += d_archive_type;
  
  ParticleSimulationManagerFactory tmp_factory( d_model,
                                                d_source,
                                                d_event_handler,
                                                d_weight_windows,
                                                d_collision_forcer,
                                                d_properties,
                                                d_simulation_name,
                                                d_archive_type,
                                                d_next_history,
                                                d_rendezvous_number+1 );

  tmp_factory.saveToFile( archive_name, true );

  ++d_rendezvous_number;
}

// Print the simulation data to the desired stream
void ParticleSimulationManager::printSimulationSummary( std::ostream& os ) const
{
  d_source->printSummary( os );
  d_event_handler->printObserverSummaries( os );
}

// Log the simulation data
void ParticleSimulationManager::logSimulationSummary() const
{
  d_source->logSummary();
  d_event_handler->logObserverSummaries();
}

// Run the simulation batch
void ParticleSimulationManager::runSimulationBatch(
                                            const uint64_t batch_start_history,
                                            const uint64_t batch_end_history )
{
  // Make sure the history range is valid
  testPrecondition( batch_start_history <= batch_end_history );

  #pragma omp parallel num_threads( Utility::OpenMPProperties::getRequestedNumberOfThreads() )
  {
    // Create a bank for each thread
    ParticleBank source_bank, bank;
    
    #pragma omp for
    for( uint64_t history = batch_start_history; history < batch_end_history; ++history )
    {
      // Initialize the random number generator for this history
      Utility::RandomNumberGenerator::initialize( history );

      // Sample a particle state from the source
      try{
        d_source->sampleParticleState( source_bank, history );
      }
      CATCH_LOST_SOURCE_PARTICLE_AND_CONTINUE( source_bank );

      // Simulate the particles generated by the source first
      while( source_bank.size() > 0 )
      {
        this->simulateUnresolvedParticle( source_bank.top(), bank, true );

        source_bank.pop();
      }

      // This history only ends when the particle bank is empty
      while( bank.size() > 0 )
      {            
        this->simulateUnresolvedParticle( bank.top(), bank, false );
        
        bank.pop();
      }

      // History complete - commit all observer history contributions
      d_event_handler->commitObserverHistoryContributions();
    }
  }
}

// The name that will be used when archiving the object
const char* ParticleSimulationManager::getArchiveName() const
{
  return s_archive_name.c_str();
}

// The signal handler
/*! \details The first signal will cause the simulation to finish. The
 * second signal will cause the simulation to end without caching its state.
 */
void ParticleSimulationManager::signalHandler( int signal )
{
  static int number_of_signals_handled = 0;

  ++number_of_signals_handled;
  
  FRENSIE_LOG_WARNING( "Terminating simulation..." );
  
  if( number_of_signals_handled == 1 )
    d_end_simulation = true;

  this->exitIfRequired( number_of_signals_handled, signal );
}

// Exit if required based on signal count
void ParticleSimulationManager::exitIfRequired( const int signal_counter,
                                                const int signal ) const
{
  if( signal_counter > 1 )
    exit( signal );
}

// Check if the simulation has been ended by the user
bool ParticleSimulationManager::hasEndSimulationRequestBeenMade() const
{
  return d_end_simulation;
}
  
} // end MonteCarlo namespace

BOOST_SERIALIZATION_CLASS_EXPORT_IMPLEMENT( MonteCarlo::ParticleSimulationManager );
EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo::ParticleSimulationManager );

//---------------------------------------------------------------------------//
// end MonteCarlo_ParticleSimulationManager.cpp
//---------------------------------------------------------------------------//
