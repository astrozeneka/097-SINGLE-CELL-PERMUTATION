#!/bin/bash

# Next.js Docker Manager Script

case "$1" in
  start)
    echo "Starting Next.js app..."
    docker-compose up -d
    echo "App started! Access it at http://localhost:3001"
    ;;
  
  stop)
    echo "Stopping Next.js app..."
    docker-compose down
    echo "App stopped."
    ;;
  
  restart)
    echo "Restarting Next.js app..."
    docker-compose restart
    echo "App restarted."
    ;;
  
  status)
    echo "=== Container Status ==="
    docker-compose ps
    echo ""
    echo "=== Health Check ==="
    docker inspect --format='{{.State.Health.Status}}' nextjs-app 2>/dev/null || echo "No health info available"
    ;;
  
  logs)
    echo "=== Recent Logs ==="
    docker-compose logs --tail=50 -f
    ;;
  
  rebuild)
    echo "Rebuilding and restarting Next.js app..."
    docker-compose down
    docker-compose build --no-cache
    docker-compose up -d
    echo "App rebuilt and started!"
    ;;
  
  *)
    echo "Next.js Docker Manager"
    echo "Usage: $0 {start|stop|restart|status|logs|rebuild}"
    echo ""
    echo "Commands:"
    echo "  start   - Start the application"
    echo "  stop    - Stop the application"
    echo "  restart - Restart the application"
    echo "  status  - Show container status and health"
    echo "  logs    - View application logs (follow mode)"
    echo "  rebuild - Rebuild and restart the application"
    exit 1
    ;;
esac
