# Use an image with both Python and Node.js
FROM nikolaik/python-nodejs:python3.11-nodejs24

WORKDIR /app

# INSTALL ZIP
RUN apt-get update && apt-get install -y zip && rm -rf /var/lib/apt/lists/*

# Copy Python requirements first for better layer caching
COPY requirements.txt ./

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

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
