FROM node:24.9.0-alpine

WORKDIR /app

# Copy package files
COPY package*.json ./

# Install all dependencies (including devDependencies needed for build)
RUN npm ci

# Copy application code
COPY . .

# Build the Next.js app
RUN npm run build

# Remove devDependencies to reduce image size
RUN npm prune --production

# Expose port 3001
EXPOSE 3001

# Set environment variable for port
ENV PORT=3001

# Start the application
CMD ["npm", "run", "start"]
